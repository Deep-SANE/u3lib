---
name: u9recoupling
description: Compute four-irrep SU(3) 9-(lambda,mu) recoupling with all six multiplicity axes and correct result lifetime (C++).
---

# u9recoupling (C++)

Read [shared runtime instructions](../runtime.md) for the execution setup and build flags.

Routes: `(A B)AB (C D)CD -> T` and `(A C)AC (B D)BD -> T`. Compute all six multiplicities with `u3::Rho` (all must be positive), then call `u3::u3u9db().get_U3U9(...)`.
Accessor: `block(iAC, iBD, iLeft, iAB, iCD, iRight)`.
* **Important Axis Ordering:** The 3rd axis is the left route `(AB)(CD) -> T`, and the 6th is the right route `(AC)(BD) -> T`. Keep them distinct, especially when their sizes differ.

**Result lifetime.** With `U3LIB_ENABLE_U3U9_CACHE` enabled (the default), `get_U3U9` returns a span-backed view into a scratch buffer that the next U9 database call overwrites. The cache-disabled build returns owning storage. To preserve a block across subsequent database calls (batching), copy its values into an owning `u3::U3U9<u3::coeffs_t>`.

```cpp
#include <u3lib.h>
#include <cassert>
#include <cmath>
#include <utility>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  // Four adjoint irreps (1,1) coupling to total scalar (0,0) via (1,1):
  // Demonstrates non-trivial outer multiplicities (rho = 2 on four axes).
  u3::irrep A(1, 1), B(1, 1), C(1, 1), D(1, 1);
  u3::irrep AB(1, 1), CD(1, 1), AC(1, 1), BD(1, 1), T(0, 0);

  int rAB = u3::Rho(A, B, AB);       // 2
  int rCD = u3::Rho(C, D, CD);       // 2
  int rAC = u3::Rho(A, C, AC);       // 2
  int rBD = u3::Rho(B, D, BD);       // 2
  int rLeft = u3::Rho(AB, CD, T);    // 1
  int rRight = u3::Rho(AC, BD, T);   // 1
  assert(rAB == 2 && rCD == 2 && rAC == 2 && rBD == 2);

  // Getter arguments: (A, B, C, D, AB, CD, AC, BD, T, rAC, rBD, rLeft, rAB, rCD, rRight)
  const auto block = u3::u3u9db().get_U3U9(A, B, C, D, AB, CD, AC, BD, T,
                                           rAC, rBD, rLeft, rAB, rCD, rRight);
  assert(block.size() == 16);

  // Verified coefficients: block(iAC, iBD, iLeft, iAB, iCD, iRight)
  assert(std::abs(block(0, 0, 0, 0, 0, 0) - (-1.0 / 3.0)) < 1e-10);
  assert(std::abs(block(0, 1, 0, 1, 0, 0) - (2.0 / 3.0)) < 1e-10);
  assert(std::abs(block(1, 1, 0, 1, 1, 0) - 0.2) < 1e-10);

  // Preserve coefficients in an owning container before any subsequent DB calls:
  u3::coeffs_t values(block.data().begin(), block.data().end());
  u3::U3U9<u3::coeffs_t> saved(std::move(values), rAC, rBD, rLeft, rAB, rCD, rRight);
  assert(std::abs(saved(1, 1, 0, 1, 1, 0) - 0.2) < 1e-10);

  u3::finalize();
}
```

For a complete recoupling matrix:

1. Discover left-route pairs `(AB, CD)` using `u3::decomposition(A, B)` and
   `u3::decomposition(C, D)`; retain pairs with `u3::Rho(AB, CD, T) > 0`.
   Enumerate every `(iAB, iCD, iLeft)` for each retained pair.
2. Discover right-route pairs `(AC, BD)` using `u3::decomposition(A, C)` and
   `u3::decomposition(B, D)`; retain pairs with `u3::Rho(AC, BD, T) > 0`.
   Enumerate every `(iAC, iBD, iRight)` for each retained pair.
3. Label matrix rows by left-route channels and columns by right-route channels.
   Fill each entry from `block(iAC, iBD, iLeft, iAB, iCD, iRight)` for the
   corresponding four intermediate irreps. Preserve the library multiplicity basis.
4. Check row and column orthogonality only after assembling both complete channel
   sets. A block for fixed intermediate irreps is generally not orthogonal by itself.
See [test/recoupling_cases.h](../../../test/recoupling_cases.h) for full matrix checks.

### Parallel batch (OpenMP)

```cpp
#include <omp.h>
#include <u3lib.h>
#include <iostream>
#include <utility>
#include <vector>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  struct U9Task {
    u3::irrep A, B, C, D, AB, CD, AC, BD, T;
    int rAC, rBD, rLeft, rAB, rCD, rRight;
  };
  u3::irrep adj(1, 1), sc(0, 0);
  std::vector<U9Task> batch = {
    {adj, adj, adj, adj, sc, sc, sc, sc, sc, 1, 1, 1, 1, 1, 1},
    {adj, adj, adj, adj, adj, adj, adj, adj, sc, 2, 2, 1, 2, 2, 1}
  };
  std::vector<double> results(batch.size());

#pragma omp parallel
  {
    u3::initialize_thread();
#pragma omp for
    for (int i = 0; i < static_cast<int>(batch.size()); ++i) {
      const auto& t = batch[i];
      const auto view = u3::u3u9db().get_U3U9(t.A, t.B, t.C, t.D,
                                              t.AB, t.CD, t.AC, t.BD, t.T,
                                              t.rAC, t.rBD, t.rLeft,
                                              t.rAB, t.rCD, t.rRight);
      // Copy span-backed view into owning container before any subsequent DB calls
      u3::coeffs_t values(view.data().begin(), view.data().end());
      u3::U3U9<u3::coeffs_t> saved(std::move(values), t.rAC, t.rBD, t.rLeft,
                                   t.rAB, t.rCD, t.rRight);
#pragma omp critical
      results[i] = saved(0, 0, 0, 0, 0, 0);
    }
    u3::finalize_thread();
  }

  for (size_t i = 0; i < results.size(); ++i)
    std::cout << "[" << i << "] = " << results[i] << "\n";

  u3::finalize();
}
```
