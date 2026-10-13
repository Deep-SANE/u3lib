---
name: decomposition
description: Enumerate the full SU(3)/U(3) Kronecker product series A x B, listing every output irrep with its outer multiplicity, using u3lib (C++).
---

# Tensor-product (Clebsch-Gordan) decomposition

Read [shared runtime instructions](../runtime.md). Use the free functions in
[U3decomposition.h](../../../include/U3decomposition.h):

API signatures (reference only):

```text
std::pair<std::vector<su3::irrep>, std::vector<int>>
su3::decomposition(const su3::irrep& A, const su3::irrep& B);   // SU(3)

std::pair<std::vector<u3::irrep>, std::vector<int>>
u3::decomposition(const u3::irrep& A, const u3::irrep& B);      // U(3)
```

`u3::decomposition` returns normalized U(3) representatives of the SU(3)
product content; the returned U(3) labels do not preserve the total U(1) charge.

## Labels, ordering and conventions

- `su3::decomposition` takes/returns SU(3) `(lambda,mu)` irreps. `u3::decomposition`
  takes U(3) irreps and returns U(3) irreps `normalize()`d (`n33=0`); it also
  reorders its two inputs by dimension internally, which does not change the set.
  Neither output list retains determinant shifts. To report full U(3) outputs,
  keep `N=A.C1()+B.C1()` and restore each output with
  `s=(N-C.lam()-2*C.mu())/3`, giving `[C.lam()+C.mu()+s,C.mu()+s,s]`.
  Check that s is a nonnegative integer and fits the library label types.
- The returned order is the generation order, **not** sorted by dimension or
  labels. Present every `(C, rhomax)` pair with explicit labels; do not assume a
  position.
- `rhomax` is a count (outer multiplicity); valid rho indices are `0..rhomax-1`
  when you later request coupling coefficients. A pair appears once with its full
  multiplicity, not repeated `rhomax` times.
- Do not cap `lambda,mu` at 5: that bound is only the regression-test range, not
  an API limit. Large products can be sizable; inspect the label types before
  unusually large inputs.

## Executable driver and checks

Drop this in `src/` (or compile against the library's configured flags, per the
runtime) and run it. The dimension sum rule is the decisive check.

```cpp
#include <u3lib.h>
#include <cassert>
#include <vector>
#include <U3decomposition.h>
#include <cstddef>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);
  su3::irrep A(1, 1), B(1, 1);
  auto [Cs, rhos] = su3::decomposition(A, B);

  // (1,1) x (1,1) = (3,0) + (2,2) + 2(1,1) + (0,3) + (0,0).
  long total = 0;
  int rho_11 = 0;
  for (std::size_t i = 0; i < Cs.size(); ++i)
  {
    total += static_cast<long>(rhos[i]) * Cs[i].dim();
    if (Cs[i] == su3::irrep(1, 1)) rho_11 = rhos[i];
    assert(rhos[i] == su3::Rho(A, B, Cs[i]));   // matches the outer-mult skill
  }
  assert(rho_11 == 2);                          // (1,1) appears with multiplicity 2
  assert(total == A.dim() * B.dim());           // sum rhomax*dim(C) = 8*8 = 64
  u3::finalize();
}
```

## Presentation and checks

Return the list of `(C, rhomax)` with labels, and state whether labels are SU(3)
`(lambda,mu)` or U(3) labels. Always report `sum_C rhomax*dim(C) == dim(A)*dim(B)`;
a mismatch means the enumeration is incomplete or an input was invalid. The set
is symmetric under `A<->B` and the conjugate series is the irrep-conjugate set.
Do not drop a `rhomax>1` entry to one copy, and do not silently restrict the
output labels.
