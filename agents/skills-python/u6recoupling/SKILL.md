---
name: u6recoupling
description: Compute U-type three-irrep recoupling blocks or complete multiplicity-resolved matrices (Python).
---

# u6recoupling

Read [shared runtime instructions](../runtime.md) before execution.

Routes: `(A B)AB C -> T` and `A (B C)BC -> T`.
Construct each irrep from the user's labels and calculate all four counts with `u3.Rho`; proceed only if all are positive:

```python
import u3lib as u3

u3.initialize()
try:
    # Example with non-trivial outer multiplicities (rho = 2 on all 4 channels):
    # Three adjoints (1,1) coupling to (1,1) through intermediate (1,1):
    A = B = C = T = u3.irrep(1, 1)
    AB = BC = u3.irrep(1, 1)

    rAB   = u3.Rho(A, B, AB)     # 2
    rAB_C = u3.Rho(AB, C, T)     # 2
    rBC   = u3.Rho(B, C, BC)     # 2
    rA_BC = u3.Rho(A, BC, T)     # 2
    assert rAB == 2 and rAB_C == 2 and rBC == 2 and rA_BC == 2

    # Note: T precedes C in get_U3U6!
    block = u3.u3u6db().get_U3U6(A, B, T, C, AB, rAB, rAB_C, BC, rBC, rA_BC)
    assert block.size() == 16
    assert abs(block(0, 0, 0, 0) - (2.0 / 9.0)) < 1e-10
    assert abs(block(0, 0, 0, 1) - 0.496903995) < 1e-6
    assert abs(block(1, 1, 1, 1) - (34.0 / 45.0)) < 1e-10

    # Dynamic channel discovery pattern for complete matrix assembly:
    irreps_AB, mults_AB = u3.decomposition(A, B)
    left_channels = [
        (int(X.lam()), int(X.mu()), e1, e2)
        for X, m in zip(irreps_AB, mults_AB)
        for e1 in range(m)
        for e2 in range(u3.Rho(X, C, T))
    ]
    assert len(left_channels) > 0
finally:
    u3.finalize()
```

Accessor: `block(iAB, iAB_C, iBC, iA_BC)`, each index bounded by its count.
* **Important Argument Order:** `T` precedes `C` in `get_U3U6(A, B, T, C, ...)`.

For a complete transformation matrix:
1. Enumerate left route channels `(AB, rho_AB, rho_AB_C)` using `u3.decomposition(A, B)` and `u3.Rho(AB, C, T)`.
2. Enumerate right route channels `(BC, rho_BC, rho_A_BC)` using `u3.decomposition(B, C)` and `u3.Rho(A, BC, T)`.
3. Sort channels lexicographically and verify row/column orthogonality.
