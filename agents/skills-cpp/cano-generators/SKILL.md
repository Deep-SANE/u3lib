---
name: cano-generators
description: Compute canonical Gelfand-state generator Eij matrix elements and all resulting states using u3lib (C++).
---

# Canonical generator action

Read [shared runtime instructions](../runtime.md). Use the methods on
`u3::gelfand` in [gelfand.h](../../../include/internal/gelfand/gelfand.h).
These are full canonical matrix elements, not SO(3)-reduced generator matrix elements.

For source g=(n13,n23,n33,n12,n22,n11), report
`Eij |g> = sum_h coefficient(h,g) |h>`; coefficient(h,g)=<h|Eij|g>.
The U(3) labels `(n13,n23,n33)` remain unchanged. Validate all six integers and interlacing before
constructing the source. Do not normalize away a determinant shift: it affects
E11,E22,E33.

## Target states and source-state methods

Each delta below is applied to (n12,n22,n11). Evaluate the named method on the
**source** state, not the target. Suffixes identify separate target branches.

| Generator | Method on g | Delta |
|---|---|---|
| E11 | e11() | (0,0,0) |
| E22 | e22() | (0,0,0) |
| E33 | e33() | (0,0,0) |
| E12 | e12() | (0,0,+1) |
| E21 | e21() | (0,0,-1) |
| E23 | e23_12() | (+1,0,0) |
| E23 | e23_22() | (0,+1,0) |
| E32 | e32_12() | (-1,0,0) |
| E32 | e32_22() | (0,-1,0) |
| E13 | e13_12() | (+1,0,+1) |
| E13 | e13_22() | (0,+1,+1) |
| E31 | e31_12() | (-1,0,-1) |
| E31 | e31_22() | (0,-1,-1) |

Reject invalid target patterns **before evaluating the method**. Some branch
formulas contain a denominator n12-n22: calling a nonexistent branch on a
boundary state can produce NaN. Do not evaluate everything and then filter NaN.
Keep the signs returned by the library, including the negative e13_22/e31_22
branches. An action can yield zero, one, or two nonzero terms.

## Executable action example

The Gelfand constructor stores unsigned 8-bit labels; validate them in the
range 0..255 before construction. This example checks the two allowed E13
branches of one source. To handle arbitrary sources, use the table above to
validate each target before calling its source-state method.

```cpp
#include <u3lib.h>
#include <cassert>
#include <cmath>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);
  u3::gelfand g(2, 1, 0, 1, 0, 0);                        // source state
  assert(std::abs(g.e13_12() - std::sqrt(0.5)) < 1e-12);  // E13 -> (2,1,0,2,0,1)
  assert(std::abs(g.e13_22() + std::sqrt(1.5)) < 1e-12);  // E13 -> (2,1,0,1,1,1), negative
  u3::finalize();
}
```

## Presentation and checks

Return each target pattern with its signed coefficient; an empty list means
zero action, not a missing result. For a requested single <h|Eij|g>, extract
the term matching h (zero if absent). For a superposition, apply each source action,
multiply by its amplitude and accumulate coincident targets without conjugating
amplitudes; generator action is linear.

For matrices, enumerate `u3::gelfandBasis(irrep)` and use target rows, source
columns. Check Eij transpose = Eji in this real basis, sum_i Eii = C1 I, and
[Eij,Ekl] = delta(j,k) Eil - delta(i,l) Ekj. Small scalar, fundamental, adjoint,
and determinant-shifted cases exercise boundary states and both branches.
Do not hide a failed identity by changing signs or discarding nonfinite values.
