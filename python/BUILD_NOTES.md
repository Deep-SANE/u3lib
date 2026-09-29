# u3lib Python bindings — build notes

`pyu3lib.cpp` builds into `u3lib.so` (a pybind11 extension). Tested on
macOS/arm64 with the GNU toolchain (Homebrew gcc-16) + conda Python 3.14.

## Building
No dependency paths are hardcoded in the Python build — they come from the
**qcx-config environment**, exactly like the C++ build. Source your env first,
then build from the u3lib root:

```
source "$QCX_ROOT/qcx-config/env.d/<you>.env.sh"
make python          # or: make -C python   (./build.sh is a wrapper for this)
```

Note: plain `make` builds only `src/` (the C++ programs); the Python module is
the separate `make python` target. `make clean` cleans it.

The `python/Makefile` reuses the qcx-config `make.d` fragments
(`boost.mk`, `eigen.mk`, `wigxjpf.mk`, `openmp.mk`), so Boost/Eigen/wigxjpf/
OpenMP all resolve from `$BOOST_ROOT`, `$EIGEN_DIR`, `$WIGXJPF_DIR`, etc. The
extension needs **none** of MPI/SLEPc/sns/libsam, so the Makefile resets the
app `CXXFLAGS`/`LDFLAGS` the top-level Makefile exports and re-includes only
those four fragments. (Verified: `otool -L u3lib.so` shows no SAM/SLEPc/PETSc/
MPI — only libgomp, plus static wigxjpf.)

## Toolchain (matches qcx-config)
- Compiler: **`g++-16`**, taken from `$OMPI_CXX` (the compiler `mpicxx` wraps).
  We use plain `g++` rather than `mpicxx` so the `.so` carries no MPI link.
  Override with `make CXX=...`.
- OpenMP: **ships with gcc (libgomp)** — native `-fopenmp` (via `$OMP_CXXFLAGS`/
  `$OMP_LDFLAGS`). No separate `libomp`, no `-Xpreprocessor`, no conda OpenMP.
- Python: `make` uses `python3` by default (must have pybind11); override with
  `make PYTHON=$HOME/miniforge3/bin/python`. Build and import with the **same**
  interpreter (ABI must match).
- Link: macOS uses `-undefined dynamic_lookup` (Python symbols resolved at load);
  other platforms use plain `-shared`.

## What the module exposes (names mirror the C++ API 1:1)
- `initialize()` / `finalize()` — **always call these** around use. `initialize()`
  cascades to the qcx-utils foundation (Wigner backend + SO3 CGC/Racah/Fano
  caches) **and** sets up the U(3) coefficient databases.
- `initialize_thread()` / `finalize_thread()` — per-thread, advanced/MT use.
- Irreps: `su3.irrep(lambda, mu)`; `irrep(lambda, mu)` and `irrep(n1, n2, n3)`
  (U(3), derives from `su3.irrep`). `str()` matches C++ `operator<<`
  (`(lambda,mu)` and `[n13,n23,n33]`).
- `gelfand(n13,n23,n33,n12,n22,n11)` — Gelfand state: labels, `check()`,
  `normalize()`, SU(2) labels, and the U(3) generator matrix elements
  (`e11`…`e31_22`).
- `Rho(...)` — coupling multiplicity (by labels or by irrep objects).
- Databases (the **only** public way to reach coefficients — there are no free
  functions): `wigU3U2db()`, `wigU3SO3db()`, `u3u6db()`, `u3z6db()`, `u3u9db()`.
- Structured results: the Definitions.h containers `U3U2`, `U3SO3`, `U3U6`,
  `U3Z6`, `U3U9` are wrapped. They give multiplicity-indexed access via
  `__call__` (e.g. `obj(kappaA, kappaB, kappaC, rho)`), flat `__getitem__`, and
  the whole buffer as numpy via `.data()`. Reached through `get_U3SO3(...)`,
  `get_U3U6(...)`, `get_U3Z6(...)`, `get_U3U9(...)`. `wigU3U2db().get(...)` and
  `wigU3SO3db().get(...)` also return flat numpy rows directly.

## Dependencies
- pybind11, Eigen (Homebrew), Boost headers (Homebrew).
- Wigner backend: **wigxjpf** (`-DUSE_WIGXJPF`, `-lwigxjpf`).
- Built **header-only** against the u3lib/qcx-utils headers (no `libu3lib`).

## Nullspace backend: QR, not sns (deliberate)
The coefficients are computed via a null-space solve. We build with the **dense
QR** backend (Eigen-only, header-only), *not* sns:
- `null_space_sns` **falls back to QR for small matrices** (`EA.size() < 10000`),
  and coupling-coefficient matrices are small — so the sns LOBPCG path is never
  taken for this workload. QR therefore gives **numerically identical
  coefficients** here.
- Keeping QR avoids dragging SLEPc/PETSc/MPI into the extension. (On the GNU
  toolchain the ABI clash that blocked sns under clang is gone — the
  `build-gnu` sns libs are libstdc++ — so sns *could* be linked, but there is no
  numerical reason to for these sizes.)

To use sns anyway (e.g. to match the production solver exactly), link
`$SNS_DIR/build-gnu/libsns_*.a` and add
`-DUSE_SNS -DNULLSPACE_SOLVER_SNS -DNULLSPACE_SPARSE`.

## Lifetime note (U3U9)
`U3U9_db::get_U3U9` returns a `std::span` view into a thread-local scratch
buffer that the next DB call overwrites. The binding **copies** it into an
owning `U3U9` before handing it to Python, so a saved `U3U9` stays valid across
later DB calls. (`U3SO3`/`U3U6`/`U3Z6` already own their buffers.)

## Using u3lib and u4lib together
u3lib and u4lib are separate extension modules, but they **share one wigxjpf
backend** in the process. `initialize()`/`finalize()` are idempotent *per module*
(guarded flags), but those flags do not coordinate the shared backend across
modules — so finalizing one library while the other is still in use tears the
shared backend out from under it and crashes the process.

Safe recipe when using both (e.g. an ENN precompute needing U(3) and U(4)
coefficients): initialize both up front, do all the work, finalize both at the
end.

```python
import u3lib as u3, u4lib as u4
u3.initialize(); u4.initialize()
#  ... all u3 and u4 coefficient work ...
u4.finalize(); u3.finalize()     # at the very end
```

Do **not** call one library's `finalize()` mid-workflow while the other is still
in use. (Within u4lib this is automatic: `u4.finalize()` already cascades to u3
and qcx_utils.) The Python `require_init()` guard is per-module, so it cannot
detect this cross-module case.

## Known caveats / TODO
- Reaching a coefficient database before `initialize()` or after `finalize()`
  raises a Python `RuntimeError` (the binding guards every DB accessor); it does
  not crash the kernel. (Exception: the cross-module teardown above.)
- **Asserts are on** (no `-DNDEBUG`): out-of-range indices or invalid labels
  `abort()` the Python process (e.g. indexing an empty/forbidden coupling).
  Before release, either build with `-DNDEBUG` or add input validation in the
  bindings so users get Python exceptions instead.
- Accuracy: spot-checked against recognizable values (e.g. `U3SO3 = -√(2/9)`,
  highest-weight `WigU3U2 = 1`, trivial `U3U6 = 1`); full numerical validation
  against reference tables is the author's to confirm.
