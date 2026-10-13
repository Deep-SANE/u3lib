# u3lib regression tests

From the u3lib root, after `make config`:

```sh
make test
```

Every invocation recompiles **all five** `test_*.cpp` executables into `bin/test/`,
then runs all of them. `make -j2 test` may compile in parallel; execution starts
only after every executable has built. No interactive input is needed.
The configured compiler, library initialization, the libwignernj Wigner backend
and the null-space
solver are inherited unchanged. Tests call the existing public library functions.
Checks remain enabled in release builds (`-UNDEBUG` and explicit runtime checks).

Each case prints `[TEST]` followed by `[PASS]` or `[FAIL]`, with its coupling or
irrep and what is checked. Reference comparisons explicitly say “against
documented values.” Coefficient values are not printed, including on failure;
a failure identifies the offending reference label or mathematical identity.
Each executable prints a summary, followed by an overall summary. Compilation
errors, exceptions, failed checks, or a crashed executable make `make test` fail.
Other executables still run after an executable fails.

## Coverage

- `test_basis.cpp`: all 36 irreps with lambda,mu <= 5; analytic dimension,
  Gelfand interlacing and ordering, determinant shifts, and SO(3) branching.
- `test_canonical.cpp`: complete allowed U(3)>U(2) reduced coefficient blocks
  for the couplings below; reconstruct full canonical CG columns for every
  output Gelfand state and check outer-multiplicity orthonormality; analytic
  scalar and conjugate-to-singlet checks.
- `test_physical.cpp`: complete allowed SO(3)-reduced blocks for the same
  couplings; orthogonality across kappa and outer multiplicities; one complete
  kappa/rho block per coupling checked by an explicit canonical-basis contraction
  at M_C=0; analytic scalar identity.
- `test_conjugation.cpp`: all 36 irreps and every allowed L,M; Gelfand phases,
  conjugation twice, single-M/all-M consistency, explicit repeated L_-lowering,
  basis orthonormality, and distinct-irrep conjugation. Self-conjugate action
  is compared with documented D matrices, including off-diagonal entries,
  and checked for M independence, symmetry, reconstruction and D squared = I.
- `test_recoupling.cpp`: scalar U6/Z6/U9 identities; complete two-channel
  U6/Z6 matrices for three fundamentals to (1,1), and U9 for four fundamentals
  to (0,2). Also tests complete multiplicity-resolved U6/Z6 matrices for
  the Python tutorial's (2,0) x (0,2) x (2,1) -> (2,1) example (5x5),
  three (1,1) irreps -> (1,1) (8x8), and (2,2) x (1,1) x (1,1) -> (2,2)
  (10x10). U9 includes four (1,1) irreps -> (0,0) (8x8) and -> (1,1)
  (32x32), plus the tutorial's scalar-leg identity. Every matrix entry is
  compared with a documented signed value; complete rows and columns are
  checked for orthogonality in the multiplicity cases. The 32x32 U9 case
  exercises nontrivial multiplicities on all six coupling axes.

All irrep labels below are SU(3) `(lambda,mu)`; normalized U(3) labels are
`[lambda+mu,mu,0]`. Coefficient couplings:

| Input A | Input B | Output C |
|---|---|---|
| (0,0) | (5,3) | (5,3) |
| (1,0) | (1,0) | (2,0), (0,1) |
| (1,0) | (0,1) | (0,0), (1,1) |
| (1,1) | (1,1) | (1,1) |
| (2,2) | (1,1) | (2,2) |
| (3,2) | (2,3) | (2,2) |
| (4,4) | (1,1) | (4,4) |
| (5,3) | (1,1) | (5,3) |
| (5,5) | (1,1) | (5,5) |
| (5,0) | (0,5) | (0,0) |

This is representative coefficient coverage, not an exhaustive Cartesian
product of all irreps. Conjugation/basis coverage is exhaustive within the bound.

## Numerical comparisons

Individual documented values use
`abs(actual - expected) <= 1e-12 + 1e-10 * abs(expected)`.
Accumulated matrix residuals use absolute Frobenius-norm tolerance `1e-9`;
the independent contraction uses absolute and relative tolerances `1e-10`.
Nonfinite results always fail. Expected signs and multiplicity ordering are
compared directly, without sign alignment or basis rotation.

RREF, Gram-Schmidt and the sign convention fix the multiplicity basis for
fixed ordering. This applies to the library's configured solvers (dense QR and
sns). Floating-point pivot thresholds can still matter for ill-conditioned
cases; a failure should be investigated, not resolved by automatically
regenerating references or loosening tolerances.

## Reference data

See [reference/README.md](reference/README.md) for labels, provenance and
analytic examples. `.csv` is not ignored by this repository, but `.dat` is used
for these whitespace-delimited, version-controlled tables.

`make test` only reads references. A separate tool can produce candidates:

```sh
make reference-tool
mkdir -p /tmp/u3-reference-candidate
./bin/test/generate_references /tmp/u3-reference-candidate
```

Generation calls library functions and runs the same mathematical identity
checks. Its output is a numerical regression baseline, not an independent
proof of the coefficients. Review candidate files and identities before
manually replacing tracked references. Do not use an existing reference
directory as the generation destination unless intentionally updating it.

To add a test, create `test_*.cpp`; the Makefile discovers it automatically.
Shared helpers live in `test_helpers.h`, while `*_cases.h` let the separate
reference tool exercise exactly the documented case selection.
