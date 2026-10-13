---
name: z6recoupling
description: Compute Z-type three-irrep SU(3) recoupling between AB-C and AC-B coupling routes (C++).
---

# z6recoupling (C++)

Read [shared runtime instructions](../runtime.md) for the execution setup and build flags.

Routes: `(A B)AB C -> T` and `(A C)AC B -> T`. Compute the four `u3::Rho` counts for these routes with `u3::Rho` (all must be positive), then call `get_U3Z6`.
* **Important Argument Order:** Note that the first two irrep arguments are **`B, A`** (the Z convention exchanges the first pair), while total `T` precedes `C`:
`get_U3Z6(B, A, T, C, AB, rAB, rAB_C, AC, rAC, rAC_B)`.
Accessor: `block(iAB, iAB_C, iAC, iAC_B)`.
The container getters `rhoBCmax()` and `rhoA_BCmax()` count the AC and AC_B
channels here, despite their U-type names. Do not add an exchange phase to the
returned coefficients unless converting to a different stated convention.

```cpp
#include <u3lib.h>
#include <cassert>
#include <cmath>
#include <vector>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  // Example with non-trivial outer multiplicities (rho = 2 on all 4 channels):
  // Three adjoints (1,1) coupling to (1,1) through intermediate (1,1):
  u3::irrep A(1, 1), B(1, 1), C(1, 1), T(1, 1);
  u3::irrep AB(1, 1), AC(1, 1);

  int rAB   = u3::Rho(A, B, AB);    // (1,1) x (1,1) -> (1,1) has rho = 2
  int rAB_C = u3::Rho(AB, C, T);    // (1,1) x (1,1) -> (1,1) has rho = 2
  int rAC   = u3::Rho(A, C, AC);    // (1,1) x (1,1) -> (1,1) has rho = 2
  int rAC_B = u3::Rho(AC, B, T);    // (1,1) x (1,1) -> (1,1) has rho = 2
  assert(rAB == 2 && rAB_C == 2 && rAC == 2 && rAC_B == 2);

  // Getter call (note the leading B, A):
  auto block = u3::u3z6db().get_U3Z6(B, A, T, C, AB, rAB, rAB_C, AC, rAC, rAC_B);
  assert(block.size() == 16);
  assert(std::abs(block(0, 0, 0, 0) - (-1.0 / 3.0)) < 1e-10);
  assert(std::abs(block(0, 1, 1, 0) - (2.0 / 3.0)) < 1e-10);
  assert(std::abs(block(1, 1, 1, 1) - 0.2) < 1e-10);

  // Right-route channel discovery pattern for Z-recoupling: (A C) -> Y, (Y B) -> T
  struct Channel { u3::irrep inter; int r1, r2; };
  std::vector<Channel> right_channels;
  auto [ac_irreps, ac_mults] = u3::decomposition(A, C);
  for (size_t i = 0; i < ac_irreps.size(); ++i) {
    int n2 = u3::Rho(ac_irreps[i], B, T);
    for (int e1 = 0; e1 < ac_mults[i]; ++e1)
      for (int e2 = 0; e2 < n2; ++e2)
        right_channels.push_back({ac_irreps[i], e1, e2});
  }
  assert(!right_channels.empty());

  u3::finalize();
}
```

For complete change-of-coupling matrix assembly:

1. Enumerate left route channels `(AB, rho_AB, rho_AB_C)` using `u3::decomposition(A, B)` and `u3::Rho(AB, C, T)`.
2. Enumerate right route channels `(AC, rho_AC, rho_AC_B)` using `u3::decomposition(A, C)` and `u3::Rho(AC, B, T)`.
3. Order channels lexicographically and check row/column orthogonality of the assembled matrix.
See [test/recoupling_cases.h](../../../test/recoupling_cases.h) for full matrix checks.

### Parallel batch (OpenMP)

```cpp
#include <omp.h>
#include <u3lib.h>
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  struct Z6Task {
    u3::irrep A, B, T, C, AB, AC;
    int rAB, rAB_C, rAC, rAC_B;
  };
  u3::irrep adj(1, 1), sc(0, 0);
  std::vector<Z6Task> batch = {
    {adj, adj, adj, adj, sc, sc, 1, 1, 1, 1},
    {adj, adj, adj, adj, adj, adj, 2, 2, 2, 2}
  };
  std::vector<double> results(batch.size());

#pragma omp parallel
  {
    u3::initialize_thread();
#pragma omp for
    for (int i = 0; i < static_cast<int>(batch.size()); ++i) {
      const auto& t = batch[i];
      // Note argument order: B precedes A, T precedes C
      auto blk = u3::u3z6db().get_U3Z6(t.B, t.A, t.T, t.C,
                                        t.AB, t.rAB, t.rAB_C,
                                        t.AC, t.rAC, t.rAC_B);
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
