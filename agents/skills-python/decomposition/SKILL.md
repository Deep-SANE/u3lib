---
name: decomposition
description: Enumerate the full SU(3)/U(3) Kronecker product series A x B, listing every output irrep with its outer multiplicity, using u3lib (Python).
---

# Tensor-product (Clebsch-Gordan) decomposition

Read [shared runtime instructions](../runtime.md). Both products are exposed in
the Python extension:

- `u3.su3.decomposition(A, B)` — SU(3); takes/returns `u3.su3.irrep(lambda,mu)`.
- `u3.decomposition(A, B)` — U(3); takes `u3.irrep(...)` and returns normalized
  U(3) irreps (`n33=0`), reordering its inputs by dimension internally.

The normalized outputs do not retain the total U(1) charge. For full U(3)
labels, keep `N=A.C1()+B.C1()` and restore each output with
`s=(N-C.lam()-2*C.mu())/3`, giving `[C.lam()+C.mu()+s,C.mu()+s,s]`.
Check that s is a nonnegative integer in the supported label range.

Each returns parallel lists: the output irreps `C` in `A x B` and their outer
multiplicities. For the multiplicity of a single fixed `C`, use
[outer-mult](../outer-mult/SKILL.md) (`u3.Rho(A,B,C)`), which equals the matching
entry.

```python
import u3lib as u3

u3.initialize()
try:
    A = B = u3.irrep(1, 1)
    outputs, rhomaxs = u3.decomposition(A, B)   # SU(3) product content; normalized U(3) outputs
    # (1,1) x (1,1) = (3,0) + (2,2) + 2(1,1) + (0,3) + (0,0)
    total = sum(r * C.dim() for C, r in zip(outputs, rhomaxs))
    assert total == A.dim() * B.dim() == 64
    assert all(r == u3.Rho(A, B, C) for C, r in zip(outputs, rhomaxs))
finally:
    u3.finalize()
```

## Ordering and checks

- The returned order is the generation order, **not** sorted by dimension or
  labels. Present every `(C, rhomax)` pair with explicit labels, and state
  whether they are SU(3) `(lambda,mu)` or U(3) labels.
- `rhomax` is a count (outer multiplicity); valid rho indices are `0..rhomax-1`.
  A pair appears once with its full multiplicity, not repeated.
- Always report `sum_C rhomax*dim(C) == dim(A)*dim(B)`; a mismatch means the
  enumeration is incomplete or an input was invalid. The set is symmetric under
  `A<->B`; the conjugate series is the irrep-conjugate set.
- Do not cap `lambda,mu` at 5 — that is only the regression-test range, not an
  API limit.
