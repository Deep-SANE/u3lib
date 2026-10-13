---
name: cano-list
description: Enumerate the canonical Gelfand-Tsetlin basis of a U(3)/SU(3) irrep with u3lib, including each state's U(2) to U(1) labels (Python).
---

# Canonical (Gelfand) basis enumeration

Read [shared runtime instructions](../runtime.md). Use `u3.gelfandBasis`, which returns the list of `u3.gelfand`
states spanning an irrep. This enumerates the basis; it does not act with a
generator (that is [cano-generators](../cano-generators/SKILL.md)) or transform
to the physical basis (that is [phys-to-cano](../phys-to-cano/SKILL.md)).

## What it returns

- `u3.gelfandBasis(irrep)` is the fixed-order list of all Gelfand patterns of
  the irrep: U(3) labels `[n13,n23,n33]` (constant), then `(n12,n22)` and `n11`.
  The enumeration order is `n12` descending `n13..n23`, then `n22` descending
  `n23..n33`, then `n11` descending `n12..n22`. Do not re-sort it: coupling and
  generator skills index states by this order.
- Each `gelfand` carries its SU(2)xU(1) subgroup labels: `epsilon()` (U(1)
  weight), `Lambdax2()` and `MLambdax2()` (twice the SU(2) angular momentum and
  its projection — doubled to stay integer), plus `check()` (betweenness) and
  the `e11()`..`e31_22()` generator matrix elements.
- Build the irrep with two SU(3) labels `u3.irrep(lambda,mu)` or three U(3) labels
  `u3.irrep(n13,n23,n33)`; validate `n13>=n23>=n33`. A nonzero determinant shift
  `n33` is retained. A uniform shift of all six labels leaves `epsilon()`,
  `Lambdax2()` and `MLambdax2()` unchanged, but changes the diagonal Eii
  eigenvalues. Keep the requested U(3) labels in the output.

## Enumerate and check

```python
import u3lib as u3

u3.initialize()
try:
    irr = u3.irrep(2, 0)
    basis = u3.gelfandBasis(irr)
    assert len(basis) == irr.dim()                 # basis size == irrep dimension
    assert all(g.check() for g in basis)           # every pattern is valid
    rows = [(str(g), g.epsilon(), g.Lambdax2(), g.MLambdax2()) for g in basis]
    # (2,0): the highest-weight state leads the list.
    assert str(basis[0]) == "2,0,0;2,0;2"
finally:
    u3.finalize()
```

## Presentation and checks

Return each state with its position index and U(2)>U(1) labels, in library order.
The strongest invariant is `len(gelfandBasis(irrep)) == irrep.dim()`; also every
`check()` is true, each `epsilon` is constant along an SU(2) `(Lambda)` multiplet,
and the `MLambdax2` values of a multiplet run in steps of 2 from `-Lambdax2` to
`+Lambdax2`. Report the fixed order as the library gives it; do not renumber to
match another enumeration convention.
