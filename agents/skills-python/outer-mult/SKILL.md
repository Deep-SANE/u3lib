---
name: outer-mult
description: Calculate outer multiplicity rho for a specified SU(3) tensor-product coupling (Python).
---

# outer-mult

Read [shared runtime instructions](../runtime.md) before execution.

Construct A,B,C as `u3.irrep(lam,mu)` and call `u3.Rho(A,B,C)`.
The six-integer overload is `u3.Rho(lamA,muA,lamB,muB,lamC,muC)`.
A zero result denotes a forbidden SU(3) output irrep, not an error for valid labels.
`Rho` uses only `(lambda,mu)` and ignores determinant shifts. For a full
U(3) product, separately require `C.C1() == A.C1() + B.C1()`; a positive
`Rho` alone does not establish that U(3) channel.

```python
import u3lib as u3

u3.initialize()
try:
    A = u3.irrep(1, 1)
    assert u3.Rho(A, A, A) == 2
finally:
    u3.finalize()
```

Return the count separately from valid zero-based rho indices. For a full
product decomposition, inspect [Python tutorial](../../../python/tutorial.ipynb)
and enumerate the complete product; do not cap intermediate labels at five
just because regression tests use that bound. Check sum_C rho*dim(C)=dim(A)*dim(B)
only when the decomposition is complete. This skill alone does not produce CGs.
