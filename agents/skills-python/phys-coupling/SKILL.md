---
name: phys-coupling
description: Compute SO(3)-reduced physical coupling coefficients with kappa and rho multiplicities (Python).
---

# phys-coupling

Read [shared runtime instructions](../runtime.md) before execution.

For each `(lam, mu, L)`, validate labels and `0 <= L <= lam + mu` before calling `u3.Kappa`.
Require positive `Kappa` for all legs, positive `Rho(A, B, C)`, and triangular condition
`abs(LA - LB) <= LC <= LA + LB`. Report forbidden channels before unsafe indexing.

```python
import u3lib as u3
import math

u3.initialize()
try:
    # (2,2) L=2 x (1,1) L=2 -> (2,2) L=2:
    # Demonstrates non-trivial inner (kappaA=2, kappaB=1, kappaC=2) and outer (rho=2) multiplicities!
    block = u3.wigU3SO3db().get_U3SO3(2, 2, 2,  1, 1, 2,  2, 2, 2)
    assert block.kappaAmax() == 2
    assert block.kappaBmax() == 1
    assert block.kappaCmax() == 2
    assert block.rhomax() == 2

    # Verified coefficients: block(kappaA, kappaB, kappaC, rho)
    assert abs(block(0, 0, 0, 0) - 0.416272707528) < 1e-9
    assert abs(block(0, 0, 0, 1) - (-0.203120282822)) < 1e-9
    assert abs(block(0, 0, 1, 0) - (-0.245145168923)) < 1e-9
    assert abs(block(0, 0, 1, 1) - (-0.502398095293)) < 1e-9
    assert abs(block(1, 0, 1, 1) - 0.245145168923) < 1e-9

    # Documented tutorial multiplicity-one check:
    assert abs(u3.wigU3SO3db().get_U3SO3(2, 0, 2, 2, 0, 2, 4, 0, 2)(0, 0, 0, 0)
               + math.sqrt(2.0 / 9.0)) < 1e-9
finally:
    u3.finalize()
```

Prefer the structured accessor `block(kappaA, kappaB, kappaC, rho)` over indexing a flat buffer.
These are reduced coefficients; full magnetic CGs additionally require contracting with the SO(3) CG factor:
`CG(LA, MA, LB, MB | LC, MC)`.

Check orthogonality over the complete required channel set; an individual
reduced block need not be orthogonal.
