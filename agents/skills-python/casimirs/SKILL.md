---
name: casimirs
description: Compute u3lib Casimir eigenvalues C1, C2 and C3 with explicit U(3) labels and SU(3) normalization (Python).
---

# Casimir eigenvalues

Read [shared runtime instructions](../runtime.md). Use the irrep methods exposed by `u3lib`; these return eigenvalues, not
matrices. See the [Python tutorial](../../../python/tutorial.ipynb).

## Which Casimir applies to which irrep

- **`C2()` and `C3()` apply to BOTH** `u3.su3.irrep(lambda,mu)` and
  `u3.irrep(...)` (the U(3) class inherits them). They are the SU(3) Casimirs and
  are **invariant under a determinant shift** (a uniform `n33` offset).
- **`C1()` is available only on `u3.irrep`**. It is the U(1) charge
  `n13+n23+n33`, not an eigenvalue determined by SU(3) labels alone.
  `u3.su3.irrep` has no `C1()` method. If only `(lambda,mu)` is supplied,
  ask for the U(3) representative or state an explicitly authorized choice.
  The two-label `u3.irrep(lambda,mu)` constructor chooses
  `[lambda+mu,mu,0]`, for which `C1=lambda+2*mu`. A specified total charge A
  implies `n33=(A-lambda-2*mu)/3`; this must be a nonnegative integer for the
  representations supported by these constructors.

## Labels and normalization

- For U(3) labels [n13,n23,n33], validate nonnegative integers and n13>=n23>=n33
  before constructing `u3.irrep(n13,n23,n33)`. For SU(3), use
  `u3.su3.irrep(lambda,mu)` (or `u3.irrep(lambda,mu)`, the normalized U(3) rep).
- The library returns

  ```text
  C1 = n13 + n23 + n33                                        (u3.irrep only)
  C2 = 2/3 * (lambda^2 + lambda*mu + mu^2 + 3*lambda + 3*mu)   (su3 and u3)
  C3 = (lambda-mu)*(lambda+2*mu+3)*(2*lambda+mu+3)/9           (su3 and u3)
  ```

  with lambda=n13-n23 and mu=n23-n33. C2 and C3 are the SU(3) Casimirs; do not
  present them as full U(3) quadratic/cubic invariants including the U(1)
  contribution. State these normalizations when comparing with another source.

## Calculate and check

```python
import u3lib as u3

u3.initialize()
try:
    # SU(3) irrep: C2 and C3 only (no C1 on u3.su3.irrep).
    s = u3.su3.irrep(2, 1)
    assert abs(s.C2() - 32/3) < 1e-12
    assert abs(s.C3() - 56/9) < 1e-12
    assert not hasattr(s, "C1")          # SU(3) labels do not fix C1

    # U(3) irrep: C1 (= total quanta) plus the same C2, C3.
    u = u3.irrep(2, 1)                    # normalized rep [3,1,0]
    assert u.C1() == 4
    assert abs(u.C2() - s.C2()) < 1e-12 and abs(u.C3() - s.C3()) < 1e-12

    # A determinant shift changes C1 by 3*s, leaves C2, C3 unchanged.
    shifted = u3.irrep(5, 3, 2)           # same (lambda,mu)=(2,1), shift n33=2
    assert shifted.C1() == 10             # 4 + 3*2
    assert shifted.C2() == u.C2() and shifted.C3() == u.C3()

    # Conjugation preserves C2 and reverses C3 (C2,C3 live on su3 and u3 alike).
    conj = u3.su3.irrep(1, 2)
    assert conj.C2() == s.C2() and conj.C3() == -s.C3()
    assert u3.su3.irrep(2, 2).C3() == 0   # self-conjugate -> C3 = 0
finally:
    u3.finalize()
```

Return a labelled row with the input irrep, whether it is SU(3)
`(lambda,mu)` or U(3) `[n13,n23,n33]`, the chosen representative when a shift was
assumed, and the requested eigenvalues; for a batch preserve input order. Do not
call `C1()` on an `u3.su3.irrep`, and do not use `u3.irrep.conjugate()` to infer a
shifted irrep's U(1) charge — that method returns a normalized irrep.
