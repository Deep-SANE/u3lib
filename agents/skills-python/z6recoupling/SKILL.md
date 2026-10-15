---
name: z6recoupling
description: Compute Z-type three-irrep recoupling between AB-C and AC-B routes (Python).
---

# z6recoupling

Read [shared runtime instructions](../runtime.md) before execution.

Routes: `(A B)AB C -> T` and `(A C)AC B -> T`.
Compute the counts for these routes with `u3.Rho`; proceed only if all are positive:

```python
import u3lib as u3

u3.initialize()
try:
    # Example with non-trivial outer multiplicities (rho = 2 on all 4 channels):
    # Three adjoints (1,1) coupling to (1,1) through intermediate (1,1):
    A = B = C = T = u3.irrep(1, 1)
    AB = AC = u3.irrep(1, 1)

    rAB   = u3.Rho(A, B, AB)     # 2
    rAB_C = u3.Rho(AB, C, T)     # 2
    rAC   = u3.Rho(A, C, AC)     # 2
    rAC_B = u3.Rho(AC, B, T)     # 2
    assert rAB == 2 and rAB_C == 2 and rAC == 2 and rAC_B == 2

    # Note: Getter starts with B, A (the Z convention exchanges the first pair):
    block = u3.u3z6db().get_U3Z6(B, A, T, C, AB, rAB, rAB_C, AC, rAC, rAC_B)
    assert block.size() == 16
    assert abs(block(0, 0, 0, 0) - (-1.0 / 3.0)) < 1e-10
    assert abs(block(0, 1, 1, 0) - (2.0 / 3.0)) < 1e-10
    assert abs(block(1, 1, 1, 1) - 0.2) < 1e-10

    # Right-route channel discovery for Z-recoupling: (A C) -> Y, (Y B) -> T
    irreps_AC, mults_AC = u3.decomposition(A, C)
    right_channels = [
        (int(Y.lam()), int(Y.mu()), e1, e2)
        for Y, m in zip(irreps_AC, mults_AC)
        for e1 in range(m)
        for e2 in range(u3.Rho(Y, B, T))
    ]
    assert len(right_channels) > 0
finally:
    u3.finalize()
```

Accessor: `block(iAB, iAB_C, iAC, iAC_B)`.
The container getters `rhoACmax()` and `rhoAC_Bmax()` count the AC and AC_B
channels, respectively. Do not add an exchange phase to the
returned coefficients unless converting to a different stated convention.
* **Important Argument Order:** The getter starts with `B, A` (`u3z6db().get_U3Z6(B, A, T, C, ...)`), despite the left route being `A x B`.

For a complete transformation matrix:

1. Enumerate left route channels `(AB, rho_AB, rho_AB_C)` using `u3.decomposition(A, B)` and `u3.Rho(AB, C, T)`.
2. Enumerate right route channels `(AC, rho_AC, rho_AC_B)` using `u3.decomposition(A, C)` and `u3.Rho(AC, B, T)`.
3. Sort channels lexicographically and verify row/column orthogonality.
