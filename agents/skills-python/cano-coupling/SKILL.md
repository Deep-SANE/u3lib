---
name: cano-coupling
description: Compute U(3) to U(2) reduced canonical coupling coefficients, including all outer multiplicities (Python).
---

# cano-coupling

Read [shared runtime instructions](../runtime.md) before execution.

Each input state is `(n13, n23, n33, n12, n22)`. Check `n13 >= n12 >= n23 >= n22 >= n33`
before calling. Verify outer multiplicity with `u3.Rho` and subgroup
selection rules before requesting coefficients.

```python
import u3lib as u3

u3.initialize()
try:
    # States a=(4,2,0, 4,2), b=(4,2,0, 4,0), c=(4,2,0, 4,2).
    # Irreps are (2,2) x (2,2) -> (2,2), which has outer multiplicity rho_max = 3.
    row = u3.wigU3U2db().get(4, 2, 0, 4, 2,  4, 2, 0, 4, 0,  4, 2, 0, 4, 2)
    assert len(row) == 3
    assert len(row) == u3.Rho(u3.irrep(2, 2), u3.irrep(2, 2), u3.irrep(2, 2))

    # Verified coefficients across the 3 outer-multiplicity rho channels:
    assert abs(row[0] - 0.774596669241) < 1e-9
    assert abs(row[1] - 0.146385010942) < 1e-9
    assert abs(row[2] - 0.0243975018237) < 1e-9
finally:
    u3.finalize()
```

The returned row is a 1-D NumPy array over `rho = 0 .. rho_max - 1`.
These are reduced coefficients; `n11` / magnetic projections are absent. If full
canonical CGs are requested, contract with the subgroup SU(2) CG using doubled
angular momenta and projections `(2*Lambda, 2*M_Lambda)`.
Highest-weight null-space, RREF, orthonormalization and leading-sign choices fix
the library rho basis. Preserve it.
