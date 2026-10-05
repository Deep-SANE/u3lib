# AGENTS.md — u3lib

Guidance for AI agents (and new contributors) working in **u3lib**: a C++
library for fast, accurate evaluation of **U(3) Wigner coupling** and **Racah
recoupling** coefficients.

- Paper: https://doi.org/10.1140/epjp/s13360-024-05501-8
- Manual (PDF, CI artifact): see `README.adoc`; LaTeX source in `doc/`.
- Sibling library **u4lib** (SU(4)/Wigner supermultiplet) mirrors this layout;
  **SAM** is the main consumer.

This repo is being prepared for public, AI-assisted release. A thin **MCP
server** wrapping the Python bindings is planned (see "Public-release direction"
at the end); keep the conventions and input-validation notes below accurate, as
that server will encode them.

## Build & test

Everything resolves from the **qcx-config environment** — no machine paths are
hardcoded. Source it first (from the QCX root), then use the make targets:

```bash
source "$QCX_ROOT/qcx-config/env.d/<you>.env.sh"
make          # build the C++ programs in src/ (-> bin/)
make test     # recompile AND run all regression tests (test/, -> bin/test/)
make python   # build the pybind11 extension -> python/u3lib.so
make clean
```

- `make test` is the regression gate: it rebuilds all five `test_*.cpp` and runs
  them, printing `[TEST]`/`[PASS]`/`[FAIL]` without dumping coefficient values.
  Checks stay on in release builds (`-UNDEBUG`). See `test/README.md` for
  coverage/tolerances and `test/reference/README.md` for saved values and
  conventions. **Run `make test` after any change to `src/` or `include/`.**
- `make python` is a **separate** target (plain `make` does not build it).
  Details, toolchain, and dependencies are in `python/BUILD_NOTES.md` — read it
  before touching the bindings. `./build.sh` is a wrapper for `make -C python`.
- Dependencies (via qcx-config `make.d`): Boost headers, Eigen, **wigxjpf**
  (Wigner backend, `-DUSE_WIGXJPF`), OpenMP. The Python extension deliberately
  carries **no** MPI/SLEPc/sns (QR null-space backend; see BUILD_NOTES).

## Layout

- `include/` — headers (header-heavy library; `u3lib.h` is the umbrella):
  - irreps & basis: `U3irrep.h`, `U3gelfand.h`, `U3basis.h`, `U3decomposition.h`,
    `U3generatorsRME.h`, `U3SO3_to_U3U2U1.h`.
  - coupling: `U3U2xU3U2.h` (U(3)⊃U(2) canonical CG, `WigU3U2`),
    `U3SO3xU3SO3.h` (physical SO(3)-reduced CG, `WigU3SO3`).
  - recoupling: `U3_U6Z6Coeffs.h` (U-coefficient `U6` and Z-coefficient `Z6`,
    both recouple **three** irreps), `U3_U9Coeffs.h` (9-(λ,μ), recouples four).
  - `Definitions.h` — the structured result containers (`U3U2`, `U3SO3`, `U3U6`,
    `U3Z6`, `U3U9`) with multiplicity-indexed access.
- `src/` — the `.cpp` implementations + `src/Makefile`.
- `test/` — regression suite (`test_*.cpp`, `*_cases.h`, `reference/`).
- `python/` — pybind11 wrapper (`pyu3lib.cpp` -> `u3lib.so`), `build.sh`,
  `BUILD_NOTES.md`, `tutorial.ipynb`.
- `bin/` — built executables (not source; do not edit).
- `doc/` — LaTeX manual (built by GitLab CI on `doc/*` changes).

## Using the library — conventions and gotchas

These apply to both the C++ API and the Python bindings (names mirror 1:1).

- **Initialize first.** Call `initialize()` before any coefficient work and
  `finalize()` at the end. `initialize()` cascades to the qcx-utils Wigner
  backend + SO(3) CGC/Racah/Fano caches **and** sets up the U(3) databases.
  `initialize_thread()`/`finalize_thread()` are the per-thread (MT) variants.
- **Coefficients come only from the DB singletons** — there are no free
  coefficient functions: `wigU3U2db()`, `wigU3SO3db()`, `u3u6db()`, `u3z6db()`,
  `u3u9db()`. They return the structured containers from `Definitions.h`.
- **Irreps:** `irrep(lambda, mu)` (SU(3) labels) or `irrep(n1, n2, n3)` (U(3)
  rows). `Rho(a, b, c)` gives the outer coupling multiplicity (by irreps or by
  labels).
- **Getter arguments are explicit and order-sensitive** — especially the
  recoupling multiplicity maxima, which sit in fixed positions interleaved with
  the irreps. Passing them in the wrong order trips an internal assertion (an
  abort, since asserts are on — see below). Full signatures are tabulated in
  "Getter signatures and multiplicity arguments" immediately below.

### Getter signatures and multiplicity arguments

`ρ(a,b,c) = Rho(a,b,c)` is the outer multiplicity of `a⊗b→c`. The coupling
(`Wig*`) getters take **no** multiplicity inputs — the output is a row/array
indexed by the multiplicity `ρ` (and `κ` for SO(3)); the recoupling getters
(`U6`/`Z6`/`U9`) take the multiplicity maxima explicitly, in the order shown.

- `wigU3U2db().get(` `n13a,n23a,n33a, n12a,n22a,` `n13b,n23b,n33b, n12b,n22b,`
  `n13c,n23c,n33c, n12c,n22c )` — U(3)⊃U(2) canonical reduced CG. Three irreps
  by Gelfand rows `[n13,n23,n33]` each with its U(2) sublabels `n12,n22`; `a⊗b→c`.
  Returns a numpy **row over ρ** (no multiplicity input).
- `wigU3SO3db().get(` `lam1,mu1,l1, lam2,mu2,l2, lam3,mu3,l3 )` — physical
  SO(3)-reduced CG, `(λ1μ1)l1 ⊗ (λ2μ2)l2 → (λ3μ3)l3`. Returns a **flat** numpy
  array indexed by `(κA,κB,κC,ρ)`. `get_U3SO3(...)` takes the same args and
  returns a structured `U3SO3` with `obj(κA,κB,κC,ρ)` access.
- `u3u6db().get_U3U6(` `A, B, ABC, C, AB,` `ρ(A,B,AB), ρ(AB,C,ABC),`
  `BC, ρ(B,C,BC), ρ(A,BC,ABC) )` — U-coefficient, recouples
  `(A·B)_{AB}·C ↔ A·(B·C)_{BC}`, both → `ABC`. The two multiplicity pairs are
  the `(AB)·C` route then the `A·(BC)` route. (Arg names: `rhoABmax, rhoAB_Cmax,
  rhoBCmax, rhoA_BCmax`.)
- `u3z6db().get_U3Z6(` `B, A, ABC, C, AB,` `ρ(A,B,AB), ρ(AB,C,ABC),`
  `AC, ρ(A,C,AC), ρ(AC,B,ABC) )` — Z-coefficient, recouples
  `(A·B)_{AB}·C ↔ (A·C)_{AC}·B`, both → `ABC`. **Note the first two irrep args
  are `B, A`** (the Z convention exchanges the first pair). Multiplicities: the
  `(AB)·C` route then the `(AC)·B` route. (Arg names: `rhoABmax, rhoAB_Cmax,
  rhoACmax, rhoAC_Bmax`.)
- `u3u9db().get_U3U9(` `A, B, C, D, AB, CD, AC, BD, ABCD,`
  `ρ(A,C,AC), ρ(B,D,BD), ρ(AB,CD,ABCD),`
  `ρ(A,B,AB), ρ(C,D,CD), ρ(AC,BD,ABCD) )` — 9-(λ,μ), recouples
  `(A·B)_{AB}·(C·D)_{CD} ↔ (A·C)_{AC}·(B·D)_{BD}`, both → `ABCD`. The **3rd**
  multiplicity is the `(AB)(CD)→ABCD` route and the **6th** is the
  `(AC)(BD)→ABCD` route — the pair most easily swapped (the swap is exactly what
  the internal `Eta/Rho(AB,CD,ABCD)` assertion catches). (Arg names: `rhoACmax,
  rhoBDmax, rhoABCDmax, rhoABmax, rhoCDmax, rhoACBDmax`.)

u4lib's getters mirror these with SU(4) irreps and `Eta` in place of `Rho`.
- **Asserts are ON (no `-DNDEBUG`).** Out-of-range indices or invalid/forbidden
  labels `abort()` the process (e.g. indexing an empty coupling, or the
  multiplicity mismatch above). In C++ this is a hard abort; in Python, DB
  access *before* `initialize()` / after `finalize()` is guarded and raises
  `RuntimeError`, but a bad label/index still aborts the interpreter. **Any
  public-facing wrapper (incl. the MCP server) must validate inputs and raise,
  not abort.**
- **`U3U9` lifetime:** `U3U9_db::get_U3U9` returns a `std::span` into a
  thread-local scratch buffer overwritten by the next DB call. The Python
  binding copies into an owning `U3U9`, so a saved Python object stays valid;
  in C++, consume or copy it before the next call. (`U3SO3`/`U3U6`/`U3Z6` own
  their buffers.)
- **Using u3lib + u4lib together:** they share one wigxjpf backend in the
  process. Initialize both up front and finalize both at the very end; do **not**
  call one library's `finalize()` while the other is still in use (it tears out
  the shared backend and crashes). See BUILD_NOTES "Using u3lib and u4lib
  together."

## Conventions & correctness

The coupling/recoupling **phase and normalization conventions** follow the paper
above; the concrete documented values and the generation procedure live in
`test/reference/` (`test/reference/README.md`). When exposing coefficients to a
new consumer (bindings, MCP server, another library), state the convention
explicitly rather than assuming — a wrong convention yields numbers that are
silently off.

Validation status: the regression suite (`make test`) checks dimensions,
Gelfand interlacing/ordering, SO(3) branching, canonical/physical CG
orthonormality across κ and outer multiplicities, conjugation identities, and
analytic scalar/conjugate checks. Python bindings are spot-checked against
recognizable values; full numerical validation against reference tables is the
author's to confirm (see BUILD_NOTES "Known caveats").

## House rules for agents

- Source qcx-config before building; never hardcode toolchain or dependency
  paths — add/adjust `make.d` fragments if something is missing.
- `bin/` and `*.dSYM/` are build output; don't commit or edit them (`.gitignore`
  covers the usual cases — check before `git add`).
- After changing `src/`/`include/`, run `make test`; after changing `python/`,
  rebuild with `make python` and import with the **same** interpreter used to
  build (ABI must match).
- Keep `README.adoc`, `BUILD_NOTES.md`, and this file in sync with the public
  API; they are the contract the release and the MCP server depend on.

## Public-release direction

Planned: a thin **MCP server** exposing u3lib (and u4lib) coupling + `U6`/`Z6`/
`U9` recoupling as tools, wrapping the existing pybind bindings. The v1 ("our
interactive use") and the public v2 share the same tool bodies; the public delta
is productization — reproducible install (incl. Linux), stable/versioned tool
schemas, docs, CI/tests — plus two things that must be right from the first
public cut: **documented, fixed coefficient conventions** and **input validation
so bad requests raise instead of aborting** (the asserts-on behavior above).
Build v1 with those two solid and it hardens into v2 without a rewrite.
