---
name: phys-to-cano
description: Expand a physical state (lambda,mu), kappa, L, M in the canonical Gelfand basis with u3lib, honoring the library's phase conventions (Python).
---

# Physical-to-canonical transformation

Read [shared runtime instructions](../runtime.md). `u3.phys_to_cano(p, M, gArr)`
expresses the physical state `|(lambda,mu) kappa L M>` in the canonical Gelfand
basis; `u3.phys_to_cano_allM(p, gArr)` does every `M = L..-L` at once. This is
the basis change [phys-coupling](../phys-coupling/SKILL.md) and
[phys-generators](../phys-generators/SKILL.md) use internally — not SU(3)
conjugation ([conjugation](../conjugation/SKILL.md)) and not a generator action
([cano-generators](../cano-generators/SKILL.md)). The highest-weight step solves
a null space.

## Inputs, output and conventions

- Use a normalized irrep (`n33=0`) and validate `0 <= L <= lambda+mu`
  before constructing `p = u3.physicalvector(irrep, L)`. Require
  `p.kappa()>0` and `|M|<=L`, then construct `gArr = u3.gelfandBasis(p.irr())`.
- `phys_to_cano(p, M, gArr)` returns `(canoBasis, T)`: `canoBasis` are indices
  into `gArr` of the Gelfand states at this `M`; **column `k` of `T` is the k-th
  SO(3) multiplicity copy** (`kappa=k+1`), row `r` the amplitude on
  `gArr[canoBasis[r]]`. `phys_to_cano_allM` returns parallel lists indexed
  `x = L - M`.
- The library phase convention is fixed (highest-weight null-space solve, then
  `L_-`, with the within-irrep `M<0` reflection relation). For a self-conjugate irrep `(n,n)`
  the multiplicity basis need not diagonalize conjugation and can mix kappa (e.g.
  `(4,4)` L=2); do not realign it. Any kappa rotation must also be applied to
  coupling coefficients.

```python
import u3lib as u3
import math

u3.initialize()
try:
    ir = u3.irrep(2, 0)
    gArr = u3.gelfandBasis(ir)
    p = u3.physicalvector(ir, 2)                  # |(2,0) L=2>, kappa_max = 1
    cb2, T2 = u3.phys_to_cano(p, 2, gArr)         # M=L=2: single highest-weight state
    assert T2.shape == (1, 1) and abs(T2[0, 0] - 1.0) < 1e-5
    cb0, T0 = u3.phys_to_cano(p, 0, gArr)         # M=0: amplitudes sqrt(1/3), sqrt(2/3)
    assert abs(T0[0, 0] - math.sqrt(1.0 / 3.0)) < 1e-5
    assert abs(T0[1, 0] - math.sqrt(2.0 / 3.0)) < 1e-5
    assert abs((T0.T @ T0)[0, 0] - 1.0) < 1e-5    # columns orthonormal: T^T T = I
    cbAll, Tall = u3.phys_to_cano_allM(p, gArr)
    assert abs(Tall[2 - 0][0, 0] - T0[0, 0]) < 1e-5
finally:
    u3.finalize()
```

Verify `T^T T = I_kappa` at each M; state this is the expansion at fixed M, not a
coupling coefficient. Do not rescale columns or drop the kappa mixing of a
self-conjugate irrep.
