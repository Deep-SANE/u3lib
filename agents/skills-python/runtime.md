# Shared execution instructions (Python)

Use these instructions with any `skills-python/` calculation skill. Read
[build notes](../../python/BUILD_NOTES.md) when setup is needed.
Shared execution conventions are given below.

## Environment and execution

Each skill's main example is a standalone calculation with the shared setup
shown below. When combining skills, combine their calculation bodies inside one
session rather than repeating their lifecycle wrappers. API signatures and mathematical formulas are reference notation, not
standalone programs.

- u3lib configures itself once with `make config` (which writes
  `config/config.generated.mk`); no environment variables are sourced. Do not
  hardcode another machine's paths. The sole Wigner backend is libwignernj.
- Use the Python interpreter matching the extension ABI. If the extension is
  missing or stale, run `make python PYTHON=/path/to/interpreter` from the repository root
  with the intended interpreter. The null-space solve uses krans (dense-QR fallback).
- Use the module path supplied by the environment and inspect `u3lib.__file__`.
  For a repository-local build, add its `python/` directory to the import path.
  Do not override an explicitly supplied runtime with a stale local extension.
- For a standalone calculation, use this lifecycle around coefficient work:

```python
import u3lib as u3
u3.initialize()
try:
    # Run the selected skill's calculation here.
    pass
finally:
    u3.finalize()
```

In an existing initialized application, respect its lifecycle instead of
finalizing a session owned by the caller. If another library sharing this
foundational state runs in the same process, initialize all of them before
joint work and finalize them only after all joint work is complete.
Use a separate process for independent calculations when ownership is unclear.

## Inputs and results

- SU(3) labels lambda,mu are nonnegative integers. U(3) labels are ordered
  nonnegative integers satisfying n13>=n23>=n33. Do not silently treat a three-label U(3) request as a two-label
  SU(3) request or remove a determinant shift without checking its meaning.
- Library integer storage is bounded: inspect the API types before unusually
  large requests. The test coverage lambda,mu<=5 is not a general API limit.
- Validate labels/indices before the call: invalid ones can abort the process.
  Do not disable assertions as an input-validation strategy.
- Multiplicity maxima are counts; valid indices are 0 through count-1. Label
  output axes explicitly, including whether the presentation is zero-based.
- Use the library's fixed basis, signs, and normalization. Do not rotate a
  multiplicity basis or align signs merely to match another convention.
- For a forbidden mathematical channel, report that fact without indexing an
  empty result. A zero entry in an allowed block is a valid coefficient.
- Return the requested numbers with labels, a reproducible script, and the check
  actually performed. Do not claim a calculation succeeded without executing it.
  Report missing dependencies if execution cannot proceed.

Use [documented references](../../test/reference/README.md) for signed regression
checks. Numerical agreement is tolerance-based, not guaranteed bitwise identity
across backends. If modifying bindings, rebuild/import with the same Python.
Never regenerate references just to make a disagreement disappear.
