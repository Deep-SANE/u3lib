---
name: phys-coupling
description: Compute SO(3)-reduced physical coupling coefficients with kappa and rho multiplicities (C++).
---

# phys-coupling (C++)

Read [shared runtime instructions](../runtime.md) before execution. The C++
examples below assume the `u3::` API.

For each `(lam, mu, L)`, validate labels and `0 <= L <= lam + mu` before calling
`u3::Kappa`. Require positive `Kappa` for all legs, positive `Rho(A, B, C)`, and
triangular condition `abs(LA - LB) <= LC <= LA + LB`. Report forbidden channels before unsafe indexing.

```cpp
#include <u3lib.h>
#include <cassert>
#include <cmath>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  // (2,2) L=2 x (1,1) L=2 -> (2,2) L=2:
  // Demonstrates non-trivial inner (kappaA=2, kappaB=1, kappaC=2) and outer (rho=2) multiplicities!
  auto block = u3::wigU3SO3db().get_U3SO3(2, 2, 2,  1, 1, 2,  2, 2, 2);
  assert(block.kappaAmax() == 2);
  assert(block.kappaBmax() == 1);
  assert(block.kappaCmax() == 2);
  assert(block.rhomax() == 2);
  assert(block.size() == 8);

  // Verified coefficients: block(kappaA, kappaB, kappaC, rho)
  assert(std::abs(block(0, 0, 0, 0) - 0.416272707528) < 1e-9);
  assert(std::abs(block(0, 0, 0, 1) - (-0.203120282822)) < 1e-9);
  assert(std::abs(block(0, 0, 1, 0) - (-0.245145168923)) < 1e-9);
  assert(std::abs(block(0, 0, 1, 1) - (-0.502398095293)) < 1e-9);
  assert(std::abs(block(1, 0, 1, 1) - 0.245145168923) < 1e-9);

  // Documented tutorial multiplicity-one check:
  assert(std::abs(u3::wigU3SO3db().get_U3SO3(2, 0, 2, 2, 0, 2, 4, 0, 2)(0, 0, 0, 0)
                  + std::sqrt(2.0 / 9.0)) < 1e-9);

  u3::finalize();
}
```

Prefer the structured accessor over guessing how the flat array reshapes.
Axes are `(kappaA, kappaB, kappaC, rho)`. Report reduced coefficients; full magnetic
CGs additionally require contracting with the SO(3) CG factor for the requested M values:
`CG(LA, MA, LB, MB | LC, MC)`.
For complete matrix orthogonality, sum over the complete channel set rather than an individual block.
See [test/physical_cases.h](../../../test/physical_cases.h).

### Parallel batch (OpenMP)

```cpp
#include <omp.h>
#include <u3lib.h>
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  struct PhysTask { int l1, m1, L1, l2, m2, L2, l3, m3, L3; };
  std::vector<PhysTask> batch = {
    {2, 0, 2,  2, 0, 2,  4, 0, 2},
    {2, 2, 2,  1, 1, 2,  2, 2, 2}
  };
  std::vector<double> results(batch.size());

#pragma omp parallel
  {
    u3::initialize_thread();
#pragma omp for
    for (int i = 0; i < static_cast<int>(batch.size()); ++i) {
      const auto& t = batch[i];
      auto blk = u3::wigU3SO3db().get_U3SO3(t.l1, t.m1, t.L1,
                                            t.l2, t.m2, t.L2,
                                            t.l3, t.m3, t.L3);
#pragma omp critical
      results[i] = blk(0, 0, 0, 0);
    }
    u3::finalize_thread();
  }

  for (size_t i = 0; i < results.size(); ++i)
    std::cout << "[" << i << "] = " << results[i] << "\n";

  u3::finalize();
}
```
