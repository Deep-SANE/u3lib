---
name: phys-generators
description: Compute SO(3)-reduced matrix elements of the physical SU(3) generators L and Q between physical basis states with u3lib (Python).
---

# Physical generator RMEs (L and Q)

Read [shared runtime instructions](../runtime.md). `u3.RME_L(bra, ket)` and
`u3.RME_Q(bra, ket)` return `<(lambda,mu) kappa' L' || T || (lambda,mu) kappa L>`
for the angular-momentum generator `T=L` and the quadrupole `T=Q`, with bra and
ket in the **same** normalized irrep (`n33=0`). Each returns a tuple `(RME, bra_kappa, ket_kappa)`:
`RME` is a `(bra.kappa() x ket.kappa())` two-dimensional NumPy array indexed by zero-based
multiplicity, and `bra_kappa`/`ket_kappa` are the zero-based kappa labels `0..kappa-1` for the
rows/columns. `RME_Q` expands into the canonical basis.

## Selection rules and conventions

- Build states with `u3.physicalBasis(irrep)` (skips zero-multiplicity L) or
  `u3.physicalvector(irrep, L)`; validate `0 <= L <= lambda+mu` before
  construction and require `kappa()>0`.
- `L` is diagonal in L: `RME_L` is zero unless `bra.L()==ket.L()`, and otherwise
  every entry of the `kappa' x kappa` block equals `sqrt(L(L+1)(2L+1))` in the current implementation.
  For orthonormal multiplicity copies, physical L should preserve kappa; flag
  the off-diagonal entries when `kappa>1` as an implementation limitation.
  Do not silently alter the returned block or claim it is a verified generator matrix.
- Define `m = bra.L() - ket.L()`. `Q` is rank 2: zero unless `|L'-L| <= 2 <= L'+L`. A `ket.L()==0` ket in the
  `m=0`/`m=1` branch divides by zero (the angular-momentum-forbidden case) —
  report it as forbidden rather than indexing the result. Keep the library's
  signs.
- The reduced-matrix-element convention carries an extra `sqrt(2L'+1)` (bra);
  state this when comparing with another source.

```python
import u3lib as u3
import math

u3.initialize()
try:
    b20 = u3.physicalBasis(u3.irrep(2, 0))        # L=2 (index 0), L=0 (index 1)
    rme, bra_kappa, ket_kappa = u3.RME_L(b20[0], b20[0])
    assert abs(rme[0, 0] - math.sqrt(30)) < 1e-4             # sqrt(L(L+1)(2L+1))
    assert abs(u3.RME_L(b20[0], b20[1])[0][0, 0]) < 1e-12    # zero off the L-diagonal
    assert abs(u3.RME_Q(b20[0], b20[0])[0][0, 0] - math.sqrt(35.0 / 3.0)) < 1e-4  # m=0
    assert abs(u3.RME_Q(b20[0], b20[1])[0][0, 0] + math.sqrt(20.0 / 3.0)) < 1e-4  # m=2
finally:
    u3.finalize()
```

Report the matrix with its row (`kappa'`) and column (`kappa`) labels; a full
`<...M'|T^k_q|...M>` is `CG(L,M;k,q|L',M') * RME / sqrt(2L'+1)` in this convention. A zero inside an allowed block is a valid matrix element.
