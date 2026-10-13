---
name: phys-list
description: Enumerate the physical SO(3)-reduced basis (kappa, L) of a U(3)/SU(3) irrep with u3lib (C++).
---

# Physical (SO(3)-reduced) basis enumeration

Read [shared runtime instructions](../runtime.md). Use `u3::physicalBasis`
from [U3basis.h](../../../include/U3basis.h), which returns the list of
`u3::physicalvector` blocks for each allowed L of an irrep. Each block
represents all kappa copies at that L. This enumerates
the SO(3) content; the per-L multiplicity count alone is [inner-mult](../inner-mult/SKILL.md),
and expanding a physical state in the Gelfand basis is [phys-to-cano](../phys-to-cano/SKILL.md).


## What it returns

- `u3::physicalBasis(irrep)` lists one `physicalvector` per angular momentum `L`
  that occurs, with `L` running `n13-n33` down to `0`. **L values with zero
  SO(3) multiplicity are omitted** — the list is not every integer L.
- Each `physicalvector` reports `L()` and `kappa()` (equal to `kappamax()`), the
  number of independent copies of that L, plus `lam()/mu()` and `n13()/n23()/n33()`.
  A physicalvector is a reduced label set, not a single magnetic state; the
  `2L+1` projections are implicit.
- Build the irrep with `u3::irrep(lambda,mu)` or `u3::irrep(n13,n23,n33)`;
  validate `n13>=n23>=n33`. Validate `0 <= L <= lambda+mu` before constructing one state with
  `u3::physicalvector(irrep, L)` and check `kappa()>0` before using it.

## Enumerate and check

```cpp
#include <u3lib.h>
#include <cassert>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);
  u3::irrep ir(4, 2);
  auto basis = u3::physicalBasis(ir);          // zero-multiplicity L already dropped
  long dsum = 0;
  for (auto& v : basis) dsum += (long)v.kappa() * (2 * v.L() + 1);
  assert(dsum == ir.dim());                    // sum kappa*(2L+1) == dim(irrep)
  assert(basis.front().L() == 6 && basis.front().kappa() == 1);  // descending L
  u3::finalize();
}
```

## Presentation and checks

Return each state with its index, `L`, and `kappa` (state kappa_max). The
dimension sum rule `sum_L kappa(L)*(2L+1) == dim(irrep)` is the primary check and
equals `u3::Kappa(irrep,L)` summed with weights. Report the library's order
(descending L) and keep omitted L values omitted — a missing L means zero
multiplicity, not an error.
