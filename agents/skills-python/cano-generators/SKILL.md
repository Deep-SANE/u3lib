---
name: cano-generators
description: Compute canonical Gelfand-state generator Eij matrix elements and all resulting states using u3lib (Python).
---

# Canonical generator action

Read [shared runtime instructions](../runtime.md). Use the methods on
`u3.gelfand` for full canonical
matrix elements, not SO(3)-reduced generator matrix elements.

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

The current GT constructor stores unsigned 8-bit labels. This helper validates
that representation before constructing states, and returns all nonzero terms
as `(target_six_tuple, matrix_element)` pairs. It uses existing library methods
rather than duplicating their square-root formulas.

```python
import u3lib as u3
import math

u3.initialize()
try:
    BRANCHES = {
        "E11": [("e11", (0,0,0))],
        "E22": [("e22", (0,0,0))],
        "E33": [("e33", (0,0,0))],
        "E12": [("e12", (0,0,1))],
        "E21": [("e21", (0,0,-1))],
        "E23": [("e23_12", (1,0,0)), ("e23_22", (0,1,0))],
        "E32": [("e32_12", (-1,0,0)), ("e32_22", (0,-1,0))],
        "E13": [("e13_12", (1,0,1)), ("e13_22", (0,1,1))],
        "E31": [("e31_12", (-1,0,-1)), ("e31_22", (0,-1,-1))],
    }

    def valid_gt(t):
        if len(t) != 6 or any(type(x) is not int or not 0 <= x <= 255 for x in t):
            return False
        a,b,c,d,e,f = t
        return a >= d >= b >= e >= c and d >= f >= e

    def generator_action(source, operator):
        source = tuple(source)
        if not valid_gt(source):
            raise ValueError("Expected a valid six-integer Gelfand pattern")
        if operator not in BRANCHES:
            raise ValueError("Expected E11, E22, E33, E12, E21, E23, E32, E13 or E31")
        g = u3.gelfand(*source)
        terms = []
        for method, delta in BRANCHES[operator]:
            target = source[:3] + tuple(source[3+i] + delta[i] for i in range(3))
            if not valid_gt(target):
                continue
            value = float(getattr(g, method)())
            if not math.isfinite(value):
                raise ArithmeticError("Nonfinite generator matrix element")
            if value != 0.0:
                terms.append((target, value))
        return terms

    source = (2,1,0,1,0,0)
    terms = generator_action(source, "E13")
    assert len(terms) == 2
    assert terms[0][0] == (2,1,0,2,0,1)
    assert abs(terms[0][1] - math.sqrt(0.5)) < 1e-12
    assert terms[1][0] == (2,1,0,1,1,1)
    assert abs(terms[1][1] + math.sqrt(1.5)) < 1e-12
    assert generator_action((0,0,0,0,0,0), "E23") == []
finally:
    u3.finalize()
```

## Presentation and checks

Return each target pattern with its signed coefficient; an empty list means
zero action, not a missing result. For a requested single <h|Eij|g>, sum terms
matching h (zero if absent). For a superposition, apply each source action,
multiply by its amplitude and accumulate coincident targets without conjugating
amplitudes; generator action is linear.

For matrices, enumerate `u3.gelfandBasis(irrep)` and use target rows, source
columns. Check Eij transpose = Eji in this real basis, sum_i Eii = C1 I, and
[Eij,Ekl] = delta(j,k) Eil - delta(i,l) Ekj. Small scalar, fundamental, adjoint,
and determinant-shifted cases exercise boundary states and both branches.
Do not hide a failed identity by changing signs or discarding nonfinite values.
