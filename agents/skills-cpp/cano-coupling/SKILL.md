---
name: cano-coupling
description: Compute U(3) to U(2) reduced canonical coupling coefficients, including all outer multiplicities (C++).
---

# cano-coupling (C++)

Read [shared runtime instructions](../runtime.md) before execution. The C++
examples below assume the `u3::` API.

Each input state is `(n13, n23, n33, n12, n22)`. Check `n13 >= n12 >= n23 >= n22 >= n33`
before any call. Check outer multiplicity `u3::Rho` and subgroup
selection using the conventions in [canonical code](../../../include/U3U2xU3U2.h);
in particular account for determinant shifts when comparing U(1) weights.
Interlacing alone does not establish an allowed coupling.

```cpp
#include <u3lib.h>
#include <cassert>
#include <cmath>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  // States a=(4,2,0, 4,2), b=(4,2,0, 4,0), c=(4,2,0, 4,2).
  // Irreps are (2,2) x (2,2) -> (2,2), which has outer multiplicity rho_max = 3.
  auto row = u3::wigU3U2db().get(4, 2, 0, 4, 2,  4, 2, 0, 4, 0,  4, 2, 0, 4, 2);
  assert(row.size() == 3);
  assert(row.size() == u3::Rho(u3::irrep(2, 2), u3::irrep(2, 2), u3::irrep(2, 2)));

  // Verified coefficients across the 3 outer-multiplicity rho channels:
  assert(std::abs(row(0) - 0.774596669241) < 1e-9);
  assert(std::abs(row(1) - 0.146385010942) < 1e-9);
  assert(std::abs(row(2) - 0.0243975018237) < 1e-9);

  u3::finalize();
}
```

Return every rho entry or validate a requested rho index against `row.size()`.
These are reduced coefficients: `n11` / magnetic projections are absent. If full
canonical CGs are requested, also contract with the subgroup SU(2) CG using the
library's doubled-angular-momentum convention; consult
[canonical tests](../../../test/canonical_cases.h) for the contraction.
Do not pass a full Gelfand six-tuple to this five-label-per-state getter.
Highest-weight null-space, RREF, orthonormalization and leading-sign choices fix
the library rho basis. Preserve it.

### Parallel batch (OpenMP)

```cpp
#include <omp.h>
#include <u3lib.h>
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  struct CanoTask { int a[5], b[5], c[5]; };
  std::vector<CanoTask> batch = {
    {{2, 0, 0, 2, 0}, {1, 0, 0, 1, 0}, {3, 0, 0, 3, 0}},
    {{4, 2, 0, 4, 2}, {4, 2, 0, 4, 0}, {4, 2, 0, 4, 2}}
  };
  std::vector<double> results(batch.size());

#pragma omp parallel
  {
    u3::initialize_thread();
#pragma omp for
    for (int i = 0; i < static_cast<int>(batch.size()); ++i) {
      const auto& t = batch[i];
      auto row = u3::wigU3U2db().get(t.a[0], t.a[1], t.a[2], t.a[3], t.a[4],
                                     t.b[0], t.b[1], t.b[2], t.b[3], t.b[4],
                                     t.c[0], t.c[1], t.c[2], t.c[3], t.c[4]);
#pragma omp critical
      results[i] = row(0);
    }
    u3::finalize_thread();
  }

  for (size_t i = 0; i < results.size(); ++i)
    std::cout << "[" << i << "] = " << results[i] << "\n";

  u3::finalize();
}
```
