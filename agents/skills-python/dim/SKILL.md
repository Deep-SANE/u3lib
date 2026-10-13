---
name: dim
description: Calculate SU(3) or U(3) irrep dimensions with u3lib (Python).
---

# dim

Read [shared runtime instructions](../runtime.md) before execution.

Validate the requested label convention before constructing the irrep.
For SU(3), use `u3.irrep(lam,mu).dim()`; for U(3), use
`u3.irrep(n13,n23,n33).dim()`. Dimension is unchanged by a determinant shift.
Check against `(lam+1)*(mu+1)*(lam+mu+2)//2`, where
lam=n13-n23 and mu=n23-n33 for a U(3) request. Return an integer.

`dim()` is available on both `u3.su3.irrep` and `u3.irrep`, and is
shift-invariant.

```python
import u3lib as u3

u3.initialize()
try:
    assert u3.irrep(4, 2).dim() == 60              # SU(3) labels
    assert u3.su3.irrep(4, 2).dim() == 60          # dim on u3.su3.irrep too
    assert u3.irrep(7, 3, 1).dim() == 60           # U(3) labels, same (4,2): shift-invariant
finally:
    u3.finalize()
```

For a requested table, preserve input order and label each row. Do not confuse
dimension with the number of distinct L values or their multiplicities.
