---
name: u6recoupling
description: Compute U-type three-irrep SU(3) recoupling blocks or complete multiplicity-resolved matrices (C++).
---

# u6recoupling (C++)

Read [shared runtime instructions](../runtime.md) for the execution setup and build flags.

Routes: `(A B)AB C -> T` and `A (B C)BC -> T`. Compute the four `u3::Rho` counts (all positive), then call `get_U3U6` — note **T precedes C** in the argument order:
`get_U3U6(A, B, T, C, AB, rAB, rAB_C, BC, rBC, rA_BC)`.
Accessor: `block(iAB, iAB_C, iBC, iA_BC)`, each index bounded by its count.

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
  u3::irrep AB(1, 1), BC(1, 1);

  int rAB   = u3::Rho(A, B, AB);    // (1,1) x (1,1) -> (1,1) has rho = 2
  int rAB_C = u3::Rho(AB, C, T);    // (1,1) x (1,1) -> (1,1) has rho = 2
  int rBC   = u3::Rho(B, C, BC);    // (1,1) x (1,1) -> (1,1) has rho = 2
  int rA_BC = u3::Rho(A, BC, T);    // (1,1) x (1,1) -> (1,1) has rho = 2
  assert(rAB == 2 && rAB_C == 2 && rBC == 2 && rA_BC == 2);

  // Getter call (note T precedes C):
  auto block = u3::u3u6db().get_U3U6(A, B, T, C, AB, rAB, rAB_C, BC, rBC, rA_BC);
  assert(block.size() == 16);
  assert(std::abs(block(0, 0, 0, 0) - (2.0 / 9.0)) < 1e-10);
  assert(std::abs(block(0, 0, 0, 1) - 0.496903995) < 1e-6);
  assert(std::abs(block(1, 1, 1, 1) - (34.0 / 45.0)) < 1e-10);

  // Dynamic channel-discovery pattern for complete matrix assembly:
  struct Channel { u3::irrep inter; int r1, r2; };
  std::vector<Channel> left_channels;
  auto [ab_irreps, ab_mults] = u3::decomposition(A, B);
  for (size_t i = 0; i < ab_irreps.size(); ++i) {
    int n2 = u3::Rho(ab_irreps[i], C, T);
    for (int e1 = 0; e1 < ab_mults[i]; ++e1)
      for (int e2 = 0; e2 < n2; ++e2)
        left_channels.push_back({ab_irreps[i], e1, e2});
  }
  assert(!left_channels.empty());

  u3::finalize();
}
```

For complete change-of-coupling matrix assembly:
1. Enumerate left route channels `(AB, rho_AB, rho_AB_C)` using `u3::decomposition(A, B)` and `u3::Rho(AB, C, T)`.
2. Enumerate right route channels `(BC, rho_BC, rho_A_BC)` using `u3::decomposition(B, C)` and `u3::Rho(A, BC, T)`.
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

  struct U6Task {
    u3::irrep A, B, T, C, AB, BC;
    int rAB, rAB_C, rBC, rA_BC;
  };
  u3::irrep adj(1, 1), sc(0, 0);
  std::vector<U6Task> batch = {
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
      auto blk = u3::u3u6db().get_U3U6(t.A, t.B, t.T, t.C,
                                        t.AB, t.rAB, t.rAB_C,
                                        t.BC, t.rBC, t.rA_BC);
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
