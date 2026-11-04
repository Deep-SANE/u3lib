#!/usr/bin/env python3
"""Freeze u3lib inputs into snapshots/{control,treatment}/workspace.

Treatment keeps u3lib/agents + the current README and manual.
Control removes the agents package and strips agent references from README and Manual.tex.
The only control/treatment difference is u3lib/agents and agent-related documentation.
"""
import hashlib
import json
import os
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH = HERE.parent
LIB = BENCH.parents[1]          # u3lib
SNAP = BENCH / 'snapshots'
COMMON = BENCH / 'tasks/common.md'


def make_environment_sh(root: Path):
    content = """#!/usr/bin/env bash
export U3LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/u3lib" && pwd)"
export PYTHONPATH="$U3LIB_DIR/python:${PYTHONPATH:-}"
export PATH="${PATH}"
"""
    (root / 'environment.sh').write_text(content)
    (root / 'environment.sh').chmod(0o755)


def main():
    if not COMMON.is_file():
        raise SystemExit(f"Missing {COMMON}")
    common_bytes = COMMON.read_bytes()

    allowed_exts = {'.h', '.hpp', '.cpp', '.cc', '.c', '.md', '.sh', '.ipynb', '.tex', '.bib', '.pdf', '.dat', '.txt', '.so', '.mk'}

    for condition in ('control', 'treatment'):
        root = SNAP / condition / 'workspace'
        if root.exists():
            shutil.rmtree(root)
        root.mkdir(parents=True)
        (root / 'common.md').write_bytes(common_bytes)
        (root / 'consumer').mkdir(parents=True, exist_ok=True)
        make_environment_sh(root)

        for p in LIB.rglob('*'):
            if not p.is_file():
                continue
            rel = p.relative_to(LIB)
            if any(x in rel.parts for x in ('.git', 'benchmarks', '__pycache__', '.ipynb_checkpoints')):
                continue
            if rel.parts[0] not in ('include', 'src', 'python', 'test', 'doc', 'agents', 'utils', 'external-libs', 'config', 'bin') and len(rel.parts) != 1:
                continue
            if condition == 'control' and rel.parts[0] == 'agents':
                continue
            if p.suffix not in allowed_exts and p.name not in ('Makefile', 'LICENSE', '.gitignore'):
                continue
            dst = root / 'u3lib' / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(p, dst)
            shutil.copymode(p, dst)

    # Clean control README.md
    control_readme = SNAP / 'control/workspace/u3lib/README.md'
    text = control_readme.read_text()
    text = text.replace('**AI-assisted computational framework for U(3)/SU(3) algebra**:',
                        'computational framework for U(3)/SU(3) algebra:')
    text = text.replace('**AI-ready computational framework for U(3)/SU(3) algebra**:',
                        'computational framework for U(3)/SU(3) algebra:')
    text = text.replace('exposed through a\nC++ API, Python bindings, and an agent-driven skill layer over all of it.',
                        'exposed through a\nC++ API and Python bindings.')
    text = text.replace('exposed through a C++ API, Python bindings, and an agent-driven skill layer over all of it.',
                        'exposed through a C++ API and Python bindings.')
    text = re.sub(r'> \*\*Note for AI Agents:\*\*.*?\n\n---\n\n', '', text, flags=re.DOTALL)
    text = re.sub(r'\*\*Installing with an AI agent\.\*\*.*?(?=- \*\*Regression Tests)', '', text, flags=re.DOTALL)
    text = re.sub(r'- `agents/`[^\n]+\n', '', text)
    if '## Documentation & AI Agent Skills' in text:
        text = text.split('## Documentation & AI Agent Skills')[0].rstrip() + '\n\n## Documentation\n\n- **Full Reference Manual:** For exhaustive mathematical derivations, phase conventions, getter signatures, and solver backends, consult [doc/Manual.pdf](doc/Manual.pdf).\n'
    elif '## Agent-assisted calculations' in text:
        text = text.split('## Agent-assisted calculations')[0]
    assert not re.search(r'\bagents?\b|\bskills?\b', text, re.I), f'Unremoved agent/skill ref in control README: {re.findall(r"agents?|skills?", text, re.I)}'
    control_readme.write_text(text.rstrip() + '\n')

    # Clean control Manual.tex
    control_manual = SNAP / 'control/workspace/u3lib/doc/Manual.tex'
    if control_manual.is_file():
        mt = control_manual.read_text()
        # Remove \section{AI-ready interface} or \section{Agent-assisted interface} until \section{Build and installation} (or \section{Conclusions} or \bibliography)
        mt = re.sub(r'\\section\{(?:Agent-assisted|AI-ready) interface\}.*?(?=\\section\{Build and installation\}|\\section\{Conclusions\}|\\section\{References\}|\\bibliographystyle)',
                    '', mt, flags=re.DOTALL)
        mt = mt.replace(
            'The C++ API is mirrored one-to-one by Python\nbindings, and a repository-local catalog of agent skills lets a code-executing\nAI assistant drive any of these capabilities from a natural-language request.\nThis document is the reference manual for the public release.',
            'The C++ API is mirrored one-to-one by Python\nbindings. This document is the reference manual for the public release.')
        mt = mt.replace(
            'The C++ API is header-only; Python bindings mirror it\none-to-one; and an agent-skill catalog exposes every capability to an AI\nassistant. The remainder documents the repository structure\n(Sec.~\\ref{sec:structure}), the classes and interface (Sec.~\\ref{sec:classes}),',
            'The C++ API is header-only; Python bindings mirror it\none-to-one. The remainder documents the repository structure\n(Sec.~\\ref{sec:structure}), the classes and interface (Sec.~\\ref{sec:classes}),')
        mt = re.sub(r'    \\item \\verb\|agents/\|[^\n]*\n', '', mt)
        mt = re.sub(r'the\s+(?:agent-assisted|AI-ready)\s+interface\s+\(Sec\.~\\ref\{sec:agents\}\),?\s*', '', mt)
        control_manual.write_text(mt)

    # Generate manifest
    manifest = {
        'version': 'u3lib-bench-v1',
        'conditions': {},
        'files': {},
        'difference': 'treatment retains u3lib/agents/ and unstripped docs; control strips agents/ and references.'
    }
    for condition in ('control', 'treatment'):
        root = SNAP / condition / 'workspace'
        hashes = {
            p.relative_to(root).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob('*')) if p.is_file()
        }
        manifest['conditions'][condition] = hashes
        manifest['files'][condition] = list(hashes)

    a, b = manifest['conditions']['control'], manifest['conditions']['treatment']
    non_doc_diffs = [k for k in a if a[k] != b.get(k) and k not in ('u3lib/README.md', 'u3lib/doc/Manual.tex')]
    assert not non_doc_diffs, f"Unexpected differences: {non_doc_diffs}"
    extra_in_treatment = set(b.keys()) - set(a.keys())
    assert all(k.startswith('u3lib/agents/') for k in extra_in_treatment), f"Non-agents extra files: {extra_in_treatment}"

    (SNAP / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(f"Successfully froze u3lib snapshots: {len(a)} control files, {len(b)} treatment files.")


if __name__ == '__main__':
    main()
