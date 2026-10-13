---
name: install
description: Build u3lib from a fresh checkout (make config -> make -> make test) before any calculation skill.
---

# install

Run this once, from the u3lib root, before any calculation skill, when the
library is not built yet. First detect the built state: if
`config/config.generated.mk` exists **and** `bin/` (for C++) or
`python/u3lib.so` (for Python) is present, u3lib is already configured — skip to
the calculation skills. u3lib needs no environment variables; `make config`
resolves every dependency and records the result in `config/config.generated.mk`.

## Prerequisites (system tools; the library does not install these)

- A C++20 compiler with OpenMP (real GCC, or Clang with libomp).
- To let `make config` download dependencies: `git`, `cmake`, and `curl` or `wget`.
- For the Python extension only: a Python interpreter with `pybind11` installed.

## 1. Configure

`make config` resolves each dependency by the chain **detect → ask → install**:
it detects a system install (standard locations and `pkg-config`), otherwise
prompts for a path, otherwise downloads the pinned upstream source and builds it
into `external-libs/` (not part of the repository).

- Interactive: `make config` — press Enter at each prompt to accept a detected
  path, or to download & build an undetected dependency.
- Non-interactive: `make config CONFIG_ARGS="--yes"` — detect, and install
  anything undetected. Give explicit paths when needed, e.g.
  `make config CONFIG_ARGS="--boost=/path --eigen=/path/eigen3 --yes"`.

Dependencies are Eigen3, Boost (headers), fmt, libwignernj, and sns. Boost's
download is large (~130 MB), so prefer a system Boost when one is available.

## 2. Build and verify

- `make` — build the C++ driver programs into `bin/`.
- `make test` — build and run the regression suite. Success prints
  `PASS: all test executables`. Treat any `FAIL` as a blocking error: report it
  and do not proceed to calculations.

## 3. Python extension (only for a Python calculation)

- `make python PYTHON=/path/to/python` — build `python/u3lib.so` with an
  interpreter that has `pybind11`. Verify the import and a known value:

  ```sh
  cd python && python3 -c "import u3lib as u; u.initialize(); print(u.irrep(4,2).dim()); u.finalize()"
  ```

  should print `60`.

## Troubleshooting

- `No configuration found. Run make config first.` — the generated config is
  missing; do step 1.
- A dependency is not found and fetching is not possible (no network): re-run
  with an explicit path, e.g. `make config CONFIG_ARGS="--boost=/path"`.
- To reset to a clean checkout and start over: `make distclean`.
- Do not hand-edit `config/config.generated.mk` or disable assertions; re-run
  `make config` to change the configuration.
