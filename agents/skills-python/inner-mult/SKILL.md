---
name: inner-mult
description: Calculate SO(3) branching multiplicity kappa for an SU(3) irrep and angular momentum L (Python).
---

# inner-mult

Read [shared runtime instructions](../runtime.md) before execution.

Here inner multiplicity means SU(3)>SO(3) branching, not a weight multiplicity.
Validate nonnegative integer L and labels. For L>lambda+mu, report zero without
calling `Kappa`: the underlying implementation asserts on that range.
Otherwise call `u3.Kappa(lam,mu,L)` (or its irrep overload).

```python
import u3lib as u3

u3.initialize()
try:
    ir = u3.irrep(4, 2)
    branches = [(L, u3.Kappa(ir, L)) for L in range(4 + 2 + 1)]
    assert u3.Kappa(ir, 2) == 2
    assert sum((2*L+1)*k for L,k in branches) == ir.dim()
finally:
    u3.finalize()
```

For a branching table, omit zero multiplicities unless requested. Each kappa
copy contains 2L+1 magnetic states. The returned count is not a kappa index.
See [Python tutorial](../../../python/tutorial.ipynb).
