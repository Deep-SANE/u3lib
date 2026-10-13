---
name: inner-mult
description: Calculate SO(3) branching multiplicity kappa for an SU(3) irrep and angular momentum L (C++).
---

# inner-mult

Read [shared runtime instructions](../runtime.md) before execution. The C++
examples below assume the `u3::` API.

Here inner multiplicity means SU(3)>SO(3) branching, not a weight multiplicity.
Validate nonnegative integer L and labels. For L>lambda+mu, report zero without
calling `Kappa`: the underlying implementation asserts on that range.
Otherwise call `su3::Kappa(lam,mu,L)` or
`u3::Kappa(u3::irrep(lam,mu),L)`. There is no three-integer `u3::Kappa` overload.

```cpp
#include <u3lib.h>
#include <cassert>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);
  u3::irrep ir(4, 2);
  assert(u3::Kappa(ir, 2) == 2);
  long dsum = 0;
  for (int L = 0; L <= 4 + 2; ++L) dsum += (long)(2 * L + 1) * u3::Kappa(ir, L);
  assert(dsum == ir.dim());
  u3::finalize();
}
```

For a branching table, omit zero multiplicities unless requested. Each kappa
copy contains 2L+1 magnetic states. The returned count is not a kappa index.
See [basis definitions](../../../include/U3basis.h).
