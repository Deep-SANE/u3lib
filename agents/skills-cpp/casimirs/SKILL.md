---
name: casimirs
description: Compute u3lib Casimir eigenvalues C1, C2 and C3 with explicit U(3) labels and SU(3) normalization (C++).
---

# Casimir eigenvalues

Read [shared runtime instructions](../runtime.md). Use the irrep methods in
[U3irrep.h](../../../include/U3irrep.h); these return eigenvalues, not matrices.

## Which Casimir applies to which irrep

- **`C2()` and `C3()` apply to BOTH** `su3::irrep(lambda,mu)` and
  `u3::irrep(...)` (the U(3) class inherits them). They are the SU(3) Casimirs and
  are **invariant under a determinant shift** (a uniform `n33` offset).
- **`C1()` is available only on `u3::irrep`**. It is the U(1) charge
  `n13+n23+n33`, not an eigenvalue determined by SU(3) labels alone.
  `su3::irrep` has no `C1()` method. If only `(lambda,mu)` is supplied,
  ask for the U(3) representative or state an explicitly authorized choice.
  The two-label `u3::irrep(lambda,mu)` constructor chooses
  `[lambda+mu,mu,0]`, for which `C1=lambda+2*mu`. A specified total charge A
  implies `n33=(A-lambda-2*mu)/3`; this must be a nonnegative integer for the
  representations supported by these constructors.

## Labels and normalization

- For U(3) labels [n13,n23,n33], validate nonnegative integers and n13>=n23>=n33
  before constructing `u3::irrep(n13,n23,n33)`. For SU(3), use
  `su3::irrep(lambda,mu)` (or `u3::irrep(lambda,mu)`, the normalized U(3) irrep).
- The library returns

  ```text
  C1 = n13 + n23 + n33                                        (u3::irrep only)
  C2 = 2/3 * (lambda^2 + lambda*mu + mu^2 + 3*lambda + 3*mu)   (su3 and u3)
  C3 = (lambda-mu)*(lambda+2*mu+3)*(2*lambda+mu+3)/9           (su3 and u3)
  ```

  with lambda=n13-n23 and mu=n23-n33. C2 and C3 are the SU(3) Casimirs; do not
  present them as full U(3) quadratic/cubic invariants including the U(1)
  contribution. State these normalizations when comparing with another source.

## Calculate and check

```cpp
#include <u3lib.h>
#include <cassert>
#include <cmath>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);
  su3::irrep s(2, 1);
  assert(std::abs(s.C2() - 32.0 / 3.0) < 1e-12);
  assert(std::abs(s.C3() - 56.0 / 9.0) < 1e-12);
  // su3::irrep has NO C1(): writing s.C1() is a compile error, not a runtime one.
  u3::irrep u(2, 1);                              // normalized rep [3,1,0]
  assert(u.C1() == 4);                            // C1 = total quanta
  u3::irrep shifted(5, 3, 2);                     // same (lambda,mu)=(2,1), n33=2
  assert(shifted.C1() == 10);                     // 4 + 3*2; C2, C3 unchanged
  assert(std::abs(su3::irrep(1, 2).C3() + s.C3()) < 1e-12);  // conjugation reverses C3
  assert(su3::irrep(2, 2).C3() == 0);             // self-conjugate -> C3 = 0
  u3::finalize();
}
```

Return a labelled row with the input irrep, whether it is SU(3)
`(lambda,mu)` or U(3) `[n13,n23,n33]`, the chosen representative when a shift was
assumed, and the requested eigenvalues; for a batch preserve input order. Do not
call `C1()` on an `su3::irrep`, and do not use `u3::irrep.conjugate()` to infer a
shifted irrep's U(1) charge — that method returns a normalized irrep.
