# Shared execution instructions (C++)

Use these instructions with any `skills-cpp/` calculation skill. The build setup
and shared conventions are given below.

## Environment and execution

Each skill's main example is a standalone calculation with the shared setup
shown below. When combining skills, combine their calculation bodies inside one
session rather than repeating their lifecycle wrappers. API signatures and mathematical formulas are reference notation, not
standalone programs.

- u3lib configures itself once with `make config` (which writes
  `config/config.generated.mk`); no environment variables are sourced. Do not
  hardcode another machine's paths. The sole Wigner backend is libwignernj.
- u3lib is header-only. Compile a driver (a program in `src/`, or a standalone
  `.cpp`) against the library's configured flags: the include paths for
  `u3lib/include`, `u3lib/include/internal`, `utils`, Eigen, Boost,
  and libwignernj, plus its link flags. The null-space solver is `krans` (its
  wrapper falls back to dense QR); do not silently change the selected solver.
- For a standalone calculation, use this lifecycle around coefficient work:

```cpp
#include <u3lib.h>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);   // Wigner backend + SO(3) caches + U(3) databases
  // Run the selected skill's calculation here.
  u3::finalize();
}
```

In an existing initialized application, respect its lifecycle instead of
finalizing a session owned by the caller. If another library sharing this
foundational state runs in the same process, initialize all of them before
joint work and finalize them only after all joint work is complete.
Use a separate process for independent calculations when ownership is unclear.

## Multi-threading (OpenMP)

Per-thread cache and scratch buffers are OpenMP `threadprivate`, so threading
is supported exclusively inside an OpenMP parallel region. Compiling with
multi-threading requires `-fopenmp` (and `-DUSE_OPENMP`) under the configured
build flags. The main thread calls `u3::initialize(&argc, &argv)` once; inside
the `#pragma omp parallel` region, every thread must call `u3::initialize_thread()`
upon entry and `u3::finalize_thread()` before exit. Distribute work across
threads with `#pragma omp for` and guard shared I/O with `#pragma omp critical`:

```cpp
#include <omp.h>
#include <u3lib.h>

u3::initialize(&argc, &argv);
#pragma omp parallel
{
  u3::initialize_thread();
  int tid = omp_get_thread_num();
  #pragma omp for
  for (int i = 0; i < n_tasks; ++i) {
    auto blk = u3::wigU3SO3db().get_U3SO3(...);
    #pragma omp critical
    std::cout << "[Thread " << tid << "] " << blk(0, 0, 0, 0) << "\n";
  }
  u3::finalize_thread();
}
u3::finalize();
```

## Inputs and results

- SU(3) labels lambda,mu are nonnegative integers. U(3) labels are ordered
  nonnegative integers satisfying n13>=n23>=n33. Do not silently treat a three-label U(3) request as a two-label
  SU(3) request or remove a determinant shift without checking its meaning.
- Library integer storage is bounded: inspect the API types before unusually
  large requests. The test coverage lambda,mu<=5 is not a general API limit.
- Validate labels/indices before the call: invalid ones abort via `assert`.
  Do not disable assertions as an input-validation strategy.
- Multiplicity maxima are counts; valid indices are 0 through count-1. Label
  output axes explicitly, including whether the presentation is zero-based.
- Use the library's fixed basis, signs, and normalization. Do not rotate a
  multiplicity basis or align signs merely to match another convention.
- For a forbidden mathematical channel, report that fact without indexing an
  empty result. A zero entry in an allowed block is a valid coefficient.
- Return the requested numbers with labels, a reproducible driver, and the check
  actually performed. Do not claim a calculation succeeded without executing it.
  Report missing dependencies if execution cannot proceed.

Use [documented references](../../test/reference/README.md) for signed regression
checks. Numerical agreement is tolerance-based, not guaranteed bitwise identity
across backends. If modifying src/ or include/, run `make test`. Never regenerate
references just to make a disagreement disappear.
