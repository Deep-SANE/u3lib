---
name: phys-generators
description: Compute SO(3)-reduced matrix elements of the physical SU(3) generators L and Q between physical basis states with u3lib (C++).
---

# Physical generator RMEs (L and Q)

Read [shared runtime instructions](../runtime.md). Use the free functions in
[U3generatorsRME.h](../../../include/U3generatorsRME.h):

API signatures (reference only):

```text
Eigen::MatrixXd u3::RME_L(u3::physicalvector& bra, u3::physicalvector& ket,
                          std::vector<int>& bra_kappa, std::vector<int>& ket_kappa);
Eigen::MatrixXd u3::RME_Q(u3::physicalvector& bra, u3::physicalvector& ket,
                          std::vector<int>& bra_kappa, std::vector<int>& ket_kappa);
```

These compute `<(lambda,mu) kappa' L' || T || (lambda,mu) kappa L>` for the
angular-momentum generator `T=L` and the quadrupole generator `T=Q`, with bra
and ket in the **same** normalized irrep (`n33=0`). The functions return the
matrix and append row/column labels to the two output vectors. Use the
[src/U3generatorsRME.cpp](../../../src/U3generatorsRME.cpp) driver as a usage
reference.

## Inputs, shape and convention

- Build states with `u3::physicalBasis(irrep)` (skips zero-multiplicity L) or
  `u3::physicalvector(irrep, L)`. Validate `0 <= L <= lambda+mu` and a positive
  `kappa()` before use; `physicalBasis` already drops empty L.
- The result is an `Eigen::MatrixXd` of shape `(bra.kappa(), ket.kappa())`,
  indexed `RME(i,j)` by **zero-based** multiplicity position. The out-parameters
  `bra_kappa`/`ket_kappa` are **appended to** with the **1-based** kappa
  labels `1..kappa()` for row/column presentation. Pass empty vectors or clear
  them before each call; the functions do not clear them for you.
- The reduced-matrix-element convention carries an extra `sqrt(2L'+1)` factor
  (bra), stated in the header. State this when comparing with another source;
  use `CG(L,M;k,q|L',M') * RME / sqrt(2L'+1)` for a full magnetic
  matrix element in this convention.

### L

`L` is a rank-1 tensor, diagonal in L: `RME_L` is zero unless `bra.L()==ket.L()`,
and otherwise every entry of the `kappa' x kappa` block equals
`sqrt(L(L+1)(2L+1))`. The library fills the **whole block** with this one value;
it does not return a kappa-diagonal matrix. This describes the current implementation, not a proof of the physical
angular-momentum matrix: for orthonormal multiplicity copies, L should preserve
kappa. Flag the off-diagonal entries when `kappa>1`; do not silently change
the library output or present it as a verified generator matrix.

### Q

`Q` is a rank-2 tensor. Define `m = bra.L() - ket.L()`; `RME_Q` is zero unless
`|m| <= 2`, and must also respect `|L'-L| <= 2 <= L'+L` (two `L=0` states give a
forbidden channel). The branches are:

| m | How the library forms it |
|---|---|
| +2 | direct, from `E12` on the canonical expansion |
| +1 | direct, from `E23_*`/`E13_*` on the canonical expansion |
|  0 | direct, diagonal in L; the returned block is symmetric |
| -1 | `-(RME_Q for +1)` transposed (phase `-1`) |
| -2 | `(RME_Q for +2)` transposed (phase `+1`) |

Domain guard: the `m=0` factor has `6*L*(2L-1)` in a denominator and the `m=1`
factor has `2L` in a denominator (`L=ket.L()`). A `ket.L()==0` ket in the `m=0`
or `m=1` branch divides by zero and yields NaN — this is exactly the
angular-momentum-forbidden case, so report it as forbidden rather than indexing
the result. Keep the library's signs (including the negative `E13_*` terms); do
not resymmetrize or flip signs to force an expected pattern.

## Executable driver and checks

Drop this in `src/` (or compile against the library's configured flags, as the
runtime describes) and run it; it reproduces the reference `src/` driver
non-interactively and self-checks against library-computed values.

```cpp
#include <u3lib.h>
#include <cassert>
#include <cmath>
#include <vector>
#include <U3basis.h>
#include <U3generatorsRME.h>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);
  auto close = [](double a, double b) { return std::abs(a - b) < 1e-4; };

  // (2,0): physical basis is L=2 (index 0) and L=0 (index 1), both kappa_max=1.
  auto b20 = u3::physicalBasis(u3::irrep(2, 0));
  std::vector<int> br, kt;

  // L is diagonal in L: <L2||L||L2> = sqrt(L(L+1)(2L+1)) = sqrt(30).
  assert(close(u3::RME_L(b20[0], b20[0], br, kt)(0, 0), std::sqrt(30.0)));
  br.clear(); kt.clear();
  // and zero off the L-diagonal.
  assert(close(u3::RME_L(b20[0], b20[1], br, kt)(0, 0), 0.0));
  br.clear(); kt.clear();

  // Q, m=0: <L2||Q||L2> = sqrt(35/3).
  assert(close(u3::RME_Q(b20[0], b20[0], br, kt)(0, 0), std::sqrt(35.0 / 3.0)));
  br.clear(); kt.clear();
  // Q, m=2: <L2||Q||L0> = -sqrt(20/3); m=-2 is its transpose (same here).
  double q2  = u3::RME_Q(b20[0], b20[1], br, kt)(0, 0); br.clear(); kt.clear();
  double qm2 = u3::RME_Q(b20[1], b20[0], br, kt)(0, 0); br.clear(); kt.clear();
  assert(close(q2, -std::sqrt(20.0 / 3.0)) && close(qm2, q2));

  // (4,2): L=2 sits at index 4 with kappa_max=2. m=0 block is 2x2 symmetric.
  auto b42 = u3::physicalBasis(u3::irrep(4, 2));
  auto q0 = u3::RME_Q(b42[4], b42[4], br, kt); br.clear(); kt.clear();
  assert(q0.rows() == 2 && q0.cols() == 2);
  assert(close(q0(0, 1), q0(1, 0)));                         // symmetric
  assert(close(q0(0, 1), -std::sqrt(80.0 / 3.0)));           // off-diagonal
  // m=-1 equals -(m=+1)^T: <L3||Q||L2> row vs <L2||Q||L3> column.
  auto qp1 = u3::RME_Q(b42[3], b42[4], br, kt); br.clear(); kt.clear(); // 1x2
  auto qm1 = u3::RME_Q(b42[4], b42[3], br, kt); br.clear(); kt.clear(); // 2x1
  for (int c = 0; c < 2; ++c) assert(close(qp1(0, c), -qm1(c, 0)));
  u3::finalize();
}
```

## Presentation and checks

Report the matrix with its row (`kappa'`) and column (`kappa`) labels and the
bra/ket `(lambda,mu) kappa_max L`, as the reference driver prints. State whether
you show reduced or full matrix elements: a full `<...M'|T^k_q|...M>` additionally
is `CG(L,M;k,q|L',M') * RME / sqrt(2L'+1)`.
For `L` verify the single value `sqrt(L(L+1)(2L+1))` and the `L'!=L` zero; for
`Q` verify `m=0` symmetry, the `m<->-m` transpose/phase relation above, and that
a forbidden channel is reported rather than returned as NaN. A zero entry inside
an allowed block is a valid matrix element. Do not regenerate references to mask
a disagreement.
