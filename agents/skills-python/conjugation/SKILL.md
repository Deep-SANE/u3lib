---
name: conjugation
description: Compute Gelfand or physical-state conjugation, including self-conjugate kappa mixing (Python).
---

# conjugation

Read [shared runtime instructions](../runtime.md) before execution.

For GT states use `g.conjugate()` and `g.conjugate_phase()` after validating
interlacing. Use normalized patterns (`n33=0`) for these SU(3) phase
relations. `conjugate()` returns a normalized pattern: applying it twice
to a shifted pattern restores its normalized version, not its original U(3)
charge. The action is `C|g> = (-1)^(n23-n12-n22-n11)|g.conjugate()>`.
For complex superpositions, conjugate amplitudes as well as mapping basis states.

For distinct irreps in the configured library phase convention:
`C|(lambda,mu) k L M> = (-1)^(lambda+mu-L+M)|(mu,lambda) k L,-M>`.
For self-conjugate irreps:
`C|(n,n) k L M> = (-1)^(L-M) sum_j D[j,k]|(n,n) j L,-M>`.
D may mix kappa; lowering alone does not fix this conjugation action.

To calculate D or apply the general operator, use `phys_to_cano` /
`phys_to_cano_allM` (see
[phys-to-cano](../phys-to-cano/SKILL.md)), apply the signed GT permutation, then
project onto the conjugate physical basis.
Use the returned canonical indices to match each conjugated Gelfand pattern
to its index in the conjugate irrep's `gelfandBasis`; do not assume the
permutation preserves enumeration order.
Do not diagonalize D unless the user explicitly requests a different basis.
Check reconstruction, norm, and double conjugation. At (4,4), L=2 the documented
D is `[[4,sqrt(33)],[sqrt(33),-4]]/7` in the current basis. Label matrix axes
as target then source; apply D to conjugated amplitudes.

## Label-level conjugation

The irrep- and Gelfand-state label maps (`conjugate()`, `conjugate_phase()`) are
available through the extension; the physical-operator `D` above is built on top of
`phys_to_cano`.

```python
import u3lib as u3

u3.initialize()
try:
    # Distinct irrep: (lambda,mu) -> (mu,lambda), with the stored phase.
    s = u3.su3.irrep(1, 2)
    assert str(s.conjugate()) == "(2,1)"
    assert s.conjugate_phase() == -1.0
    # Gelfand state: conjugation maps to a valid pattern and is an involution.
    g = u3.gelfand(2, 1, 0, 1, 0, 0)
    assert g.conjugate().check()
    assert g.conjugate().conjugate() == g         # label map is an involution
finally:
    u3.finalize()
```
