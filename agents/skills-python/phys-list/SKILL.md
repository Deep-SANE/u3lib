---
name: phys-list
description: Enumerate the physical SO(3)-reduced basis (kappa, L) of a U(3)/SU(3) irrep with u3lib (Python).
---

# Physical (SO(3)-reduced) basis enumeration

Read [shared runtime instructions](../runtime.md). Use `u3.physicalBasis`, which returns the list of
`u3.physicalvector` blocks for each allowed L of an irrep. Each block
represents all kappa copies at that L. This enumerates
the SO(3) content; the per-L multiplicity count alone is [inner-mult](../inner-mult/SKILL.md),
and expanding a physical state in the Gelfand basis is [phys-to-cano](../phys-to-cano/SKILL.md).


## What it returns

- `u3.physicalBasis(irrep)` lists one `physicalvector` per angular momentum `L`
  that occurs, with `L` running `n13-n33` down to `0`. **L values with zero
  SO(3) multiplicity are omitted** — the list is not every integer L.
- Each `physicalvector` reports `L()` and `kappa()` (equal to `kappamax()`), the
  number of independent copies of that L, plus `lam()/mu()` and `n13()/n23()/n33()`.
  A physicalvector is a reduced label set, not a single magnetic state; the
  `2L+1` projections are implicit.
- Build the irrep with `u3.irrep(lambda,mu)` or `u3.irrep(n13,n23,n33)`;
  validate `n13>=n23>=n33`. Validate `0 <= L <= lambda+mu` before constructing one state with
  `u3.physicalvector(irrep, L)` and check `kappa()>0` before using it.

## Enumerate and check

```python
import u3lib as u3

u3.initialize()
try:
    irr = u3.irrep(4, 2)
    basis = u3.physicalBasis(irr)
    content = [(v.L(), v.kappa()) for v in basis]
    assert content == [(6,1),(5,1),(4,2),(3,1),(2,2),(0,1)]   # note: no L=1
    # dimension check: sum over states of kappa*(2L+1) == dim(irrep)
    assert sum(k*(2*L+1) for L,k in content) == irr.dim()
finally:
    u3.finalize()
```

## Presentation and checks

Return each state with its index, `L`, and `kappa` (state kappa_max). The
dimension sum rule `sum_L kappa(L)*(2L+1) == dim(irrep)` is the primary check and
equals `u3.Kappa(irrep,L)` summed with weights. Report the library's order
(descending L) and keep omitted L values omitted — a missing L means zero
multiplicity, not an error.
