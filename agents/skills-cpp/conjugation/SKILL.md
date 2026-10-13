---
name: conjugation
description: Compute Gelfand or physical-state conjugation, including self-conjugate kappa mixing (C++).
---

# conjugation

Read [shared runtime instructions](../runtime.md) before execution.

For Gelfand states use `g.conjugate()` and `g.conjugate_phase()` after validating
interlacing. Use normalized patterns (`n33=0`) for these SU(3) phase
relations. `conjugate()` returns a normalized pattern: applying it twice
to a shifted pattern restores its normalized version, not its original U(3)
charge. The action is `C|g> = (-1)^(n23-n12-n22-n11)|g.conjugate()>`.
For complex superpositions, conjugate amplitudes as well as mapping basis states.

In the configured library, phase convention for distinct irreps is:
`C|(lambda,mu) k L M> = (-1)^(lambda+mu-L+M)|(mu,lambda) k L,-M>`.
For self-conjugate irreps:
`C|(n,n) k L M> = (-1)^(L-M) sum_j D[j,k]|(n,n) j L,-M>`.
D may mix kappa; lowering alone does not fix this conjugation action.

To calculate D or apply the general operator, use `phys_to_cano` /
`phys_to_cano_allM` (see
[phys-to-cano](../phys-to-cano/SKILL.md)), apply the signed Gelfand permutation,
then project onto the conjugate physical basis.
Read [transformation](../../../include/U3SO3_to_U3U2U1.h) and
[conjugation checks](../../../test/conjugation_cases.h) for sparse index mapping.
Do not diagonalize D unless the user explicitly requests a different basis.
Check reconstruction, norm, and double conjugation. At (4,4), L=2 the documented
D is `[[4,sqrt(33)],[sqrt(33),-4]]/7` in the current basis. Label matrix axes
as target then source; apply D to conjugated amplitudes.

## Label-level conjugation

The irrep- and Gelfand-state label maps (`conjugate()`, `conjugate_phase()`) are
exposed; the physical-operator `D` above is built on top of `phys_to_cano`.

```cpp
#include <u3lib.h>
#include <cassert>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);
  su3::irrep s(1, 2);
  assert(s.conjugate() == su3::irrep(2, 1));      // (lambda,mu) -> (mu,lambda)
  u3::gelfand g(2, 1, 0, 1, 0, 0);
  assert(g.conjugate().check());                  // conjugate is a valid pattern
  assert(g.conjugate().conjugate() == g);         // label map is an involution
  u3::finalize();
}
```
