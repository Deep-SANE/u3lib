---
name: outer-mult
description: Calculate outer multiplicity rho for a specified SU(3) tensor-product coupling (C++).
---

# outer-mult

Read [shared runtime instructions](../runtime.md) before execution. The C++
examples below assume the `u3::` API.

Construct A,B,C as `u3::irrep(lam,mu)` and call `u3::Rho(A,B,C)`.
The six-integer overload is `su3::Rho(lamA,muA,lamB,muB,lamC,muC)`;
`u3::Rho` takes irrep objects.
A zero result denotes a forbidden SU(3) output irrep, not an error for valid labels.
`Rho` uses only `(lambda,mu)` and ignores determinant shifts. For a full
U(3) product, separately require `C.C1() == A.C1() + B.C1()`; a positive
`Rho` alone does not establish that U(3) channel.

```cpp
#include <u3lib.h>
#include <cassert>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);
  su3::irrep A(1, 1);
  assert(su3::Rho(A, A, A) == 2);               // (1,1) appears twice in (1,1)x(1,1)
  u3::finalize();
}
```

Return the count separately from valid zero-based rho indices. For a full
product decomposition, inspect [decomposition](../../../include/U3decomposition.h)
and enumerate the complete product; do not cap intermediate labels at five
just because regression tests use that bound. Check sum_C rho*dim(C)=dim(A)*dim(B)
only when the decomposition is complete. This skill alone does not produce CGs.
