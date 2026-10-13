---
name: u9recoupling
description: Compute four-irrep 9-(lambda,mu) recoupling with all six multiplicity axes and correct result lifetime (Python).
---

# u9recoupling

Read [shared runtime instructions](../runtime.md) before execution.

Routes: `(A B)AB (C D)CD -> T` and `(A C)AC (B D)BD -> T`.
Calculate all six multiplicities before the call; all must be positive:

```python
import u3lib as u3

u3.initialize()
try:
    # Four adjoint irreps (1,1) coupling to total scalar (0,0) via (1,1):
    # Demonstrates non-trivial outer multiplicities (rho = 2 on four axes).
    A = B = C = D = u3.irrep(1, 1)
    AB = CD = AC = BD = u3.irrep(1, 1)
    T = u3.irrep(0, 0)

    rAB, rCD = u3.Rho(A, B, AB), u3.Rho(C, D, CD)
    rAC, rBD = u3.Rho(A, C, AC), u3.Rho(B, D, BD)
    rLeft, rRight = u3.Rho(AB, CD, T), u3.Rho(AC, BD, T)
    assert rAB == 2 and rCD == 2 and rAC == 2 and rBD == 2
    assert rLeft == 1 and rRight == 1

    # Getter arguments: (A, B, C, D, AB, CD, AC, BD, T, rAC, rBD, rLeft, rAB, rCD, rRight)
    block = u3.u3u9db().get_U3U9(A, B, C, D, AB, CD, AC, BD, T,
                                 rAC, rBD, rLeft, rAB, rCD, rRight)
    assert block.size() == 16

    # Verified coefficients: block(iAC, iBD, iLeft, iAB, iCD, iRight)
    assert abs(block(0, 0, 0, 0, 0, 0) - (-1.0 / 3.0)) < 1e-10
    assert abs(block(0, 1, 0, 1, 0, 0) - (2.0 / 3.0)) < 1e-10
    assert abs(block(1, 1, 0, 1, 1, 0) - 0.2) < 1e-10
finally:
    u3.finalize()
```

Accessor: `block(iAC, iBD, iLeft, iAB, iCD, iRight)`.
* **Important Axis Ordering:** The 3rd multiplicity belongs to `AB x CD -> T`; the 6th to `AC x BD -> T`. Keep these distinct, especially when their sizes differ.

`u3u9db().get_U3U9(...)` returns an owning result in Python that you can keep and index across later DB calls; the returned object already owns its coefficients.
For a complete recoupling matrix:

1. Discover left-route pairs `(AB, CD)` using `u3.decomposition(A, B)` and
   `u3.decomposition(C, D)`; retain pairs with `u3.Rho(AB, CD, T) > 0`.
   Enumerate every `(iAB, iCD, iLeft)` for each retained pair.
2. Discover right-route pairs `(AC, BD)` using `u3.decomposition(A, C)` and
   `u3.decomposition(B, D)`; retain pairs with `u3.Rho(AC, BD, T) > 0`.
   Enumerate every `(iAC, iBD, iRight)` for each retained pair.
3. Label matrix rows by left-route channels and columns by right-route channels.
   Fill each entry from `block(iAC, iBD, iLeft, iAB, iCD, iRight)` for the
   corresponding four intermediate irreps. Preserve the library multiplicity basis.
4. Check row and column orthogonality only after assembling both complete channel
   sets. A block for fixed intermediate irreps is generally not orthogonal by itself.
