#!/usr/bin/env python3
"""Autonomous coding agent loop for Gemini Flash inside a benchmark sandbox."""
import argparse
import json
import os
import subprocess
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

API_KEY = os.environ.get('GEMINI_API_KEY')
DEFAULT_MODEL = 'gemini-3-flash-preview'

TOOLS = [{
    'function_declarations': [
        {
            'name': 'run_bash',
            'description': 'Execute a shell command in the workspace directory and return its stdout, stderr, and returncode.',
            'parameters': {
                'type': 'OBJECT',
                'properties': {
                    'command': {'type': 'STRING', 'description': 'The shell command to run.'}
                },
                'required': ['command']
            }
        },
        {
            'name': 'read_file',
            'description': 'Read the entire contents of a file relative to the workspace root.',
            'parameters': {
                'type': 'OBJECT',
                'properties': {
                    'path': {'type': 'STRING', 'description': 'Relative path to the file to read.'}
                },
                'required': ['path']
            }
        },
        {
            'name': 'write_file',
            'description': 'Write or overwrite text content to a file. Creates parent directories automatically.',
            'parameters': {
                'type': 'OBJECT',
                'properties': {
                    'path': {'type': 'STRING', 'description': 'Relative path to the file to write.'},
                    'content': {'type': 'STRING', 'description': 'Text content to write.'}
                },
                'required': ['path', 'content']
            }
        },
        {
            'name': 'list_directory',
            'description': 'List filenames and directories inside a given relative path.',
            'parameters': {
                'type': 'OBJECT',
                'properties': {
                    'path': {'type': 'STRING', 'description': 'Relative path to directory (default ".").'}
                }
            }
        }
    ]
}]


def execute_tool(name: str, args: dict, cwd: Path) -> dict:
    try:
        if name == 'run_bash':
            cmd = args.get('command', '')
            proc = subprocess.run(cmd, shell=True, cwd=str(cwd), text=True,
                                  capture_output=True, timeout=180)
            return {
                'returncode': proc.returncode,
                'stdout': proc.stdout[-4000:],
                'stderr': proc.stderr[-4000:],
            }
        elif name == 'read_file':
            target = (cwd / args.get('path', '')).resolve()
            if not target.is_relative_to(cwd):
                return {'error': 'Access denied: path outside workspace'}
            if not target.is_file():
                return {'error': f'File not found: {args.get("path")}'}
            text = target.read_text(encoding='utf-8', errors='replace')
            if len(text) > 30000:
                text = text[:30000] + '\n... [truncated]'
            return {'content': text}
        elif name == 'write_file':
            target = (cwd / args.get('path', '')).resolve()
            if not target.is_relative_to(cwd):
                return {'error': 'Access denied: path outside workspace'}
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(args.get('content', ''), encoding='utf-8')
            return {'success': True, 'bytes_written': len(args.get('content', ''))}
        elif name == 'list_directory':
            rel = args.get('path', '.') or '.'
            target = (cwd / rel).resolve()
            if not target.is_relative_to(cwd):
                return {'error': 'Access denied: path outside workspace'}
            if not target.is_dir():
                return {'error': f'Directory not found: {rel}'}
            entries = []
            for child in sorted(target.iterdir()):
                entries.append(f"{child.name}/" if child.is_dir() else child.name)
            return {'entries': entries[:100]}
        else:
            return {'error': f'Unknown tool {name}'}
    except Exception as e:
        return {'error': str(e)}


def call_gemini(model: str, contents: list, system_instruction: str = None) -> dict:
    url = f"https://generativelanguage.googleapis.com/v1beta/models/{model}:generateContent?key={API_KEY}"
    payload = {
        'contents': contents,
        'tools': TOOLS,
    }
    if system_instruction:
        payload['system_instruction'] = {
            'parts': [{'text': system_instruction}]
        }

    data = json.dumps(payload).encode('utf-8')
    req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'})
    
    # Retry loop with backoff for rate limits / transient errors
    for attempt in range(8):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read().decode('utf-8'))
        except urllib.error.HTTPError as e:
            err_body = e.read().decode('utf-8')
            if e.code == 429 and attempt < 7:
                delay = 15.0
                try:
                    err_json = json.loads(err_body)
                    for d in err_json.get('error', {}).get('details', []):
                        if 'retryDelay' in d:
                            delay = float(d['retryDelay'].rstrip('s')) + 2.0
                except Exception:
                    pass
                time.sleep(delay)
                continue
            elif e.code in (500, 502, 503, 504) and attempt < 7:
                time.sleep(2 ** attempt + 2)
                continue
            raise RuntimeError(f"Gemini API HTTP {e.code}: {err_body}")
        except Exception as e:
            if attempt < 7:
                time.sleep(2 ** attempt + 1)
                continue
            raise


def run_agent(workspace: Path, prompt_text: str, model: str = DEFAULT_MODEL, max_turns: int = 30) -> dict:
    start_time = time.time()
    system_instruction = (
        "You are an expert autonomous AI software engineer. You are working in a local workspace to solve "
        "a scientific software integration task. The general contract is in common.md and the specific task is in prompt.md.\n"
        "Explore the repository and any guidance or skills available in this workspace.\n"
        "Implement the complete requested component strictly inside consumer/.\n"
        "Deliver source files, local test files, and a short README.md in consumer/ with build/use instructions.\n"
        "Use the tools provided (read_file, write_file, run_bash, list_directory) to build, test, and verify your submission before finishing."
    )

    contents = [
        {'role': 'user', 'parts': [{'text': prompt_text}]}
    ]

    total_prompt_tokens = 0
    total_candidates_tokens = 0
    total_thoughts_tokens = 0
    num_turns = 0
    is_error = False
    final_result = ""

    for turn in range(max_turns):
        num_turns += 1
        resp = call_gemini(model, contents, system_instruction)
        
        usage = resp.get('usageMetadata', {})
        total_prompt_tokens += usage.get('promptTokenCount', 0)
        total_candidates_tokens += usage.get('candidatesTokenCount', 0)
        total_thoughts_tokens += usage.get('thoughtsTokenCount', 0)

        candidates = resp.get('candidates', [])
        if not candidates:
            is_error = True
            final_result = "No candidates returned by Gemini"
            break

        cand = candidates[0]
        content = cand.get('content', {})
        parts = content.get('parts', [])
        contents.append(content)

        tool_calls = [p['functionCall'] for p in parts if 'functionCall' in p]
        text_parts = [p['text'] for p in parts if 'text' in p]

        if not tool_calls:
            final_result = "\n".join(text_parts)
            break

        # Execute tool calls
        response_parts = []
        for tc in tool_calls:
            call_id = tc.get('id')
            fn_name = tc.get('name')
            fn_args = tc.get('args', {})
            tool_res = execute_tool(fn_name, fn_args, workspace)
            
            part = {
                'functionResponse': {
                    'name': fn_name,
                    'response': {'result': tool_res}
                }
            }
            if call_id:
                part['functionResponse']['id'] = call_id
            response_parts.append(part)

        contents.append({
            'role': 'user',
            'parts': response_parts
        })

    elapsed = time.time() - start_time
    # Approximate pricing for Gemini Flash ($0.10 / 1M prompt, $0.40 / 1M candidate)
    cost = (total_prompt_tokens * 0.10 + total_candidates_tokens * 0.40) / 1_000_000.0

    return {
        'duration_ms': int(elapsed * 1000),
        'total_cost_usd': round(cost, 6),
        'usage': {
            'input_tokens': total_prompt_tokens,
            'output_tokens': total_candidates_tokens,
            'thoughts_tokens': total_thoughts_tokens,
            'total_tokens': total_prompt_tokens + total_candidates_tokens,
        },
        'num_turns': num_turns,
        'is_error': is_error,
        'result': final_result,
    }


def main():
    parser = argparse.ArgumentParser(description="Gemini autonomous benchmark agent")
    parser.add_argument('--workspace', type=Path, default=Path.cwd(), help="Workspace directory")
    parser.add_argument('--model', type=str, default=DEFAULT_MODEL, help="Gemini model name")
    parser.add_argument('--prompt', type=str,
                        default="Implement the task specified in common.md and prompt.md. Write all deliverables only in consumer/.")
    args = parser.parse_args()

    res = run_agent(args.workspace.resolve(), args.prompt, model=args.model)
    print(json.dumps(res, indent=2))
    sys.exit(0 if not res['is_error'] else 1)


if __name__ == '__main__':
    main()
