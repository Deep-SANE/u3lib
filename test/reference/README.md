# Documented coefficient and conjugation values

The `.dat` files contain `label value` rows with 17 significant digits.
Comment lines start with `#`. The labels specify quantum numbers; integers
in labels are not inferred from row position. Every documented entry must be
exercised, and missing, duplicate or nonfinite entries are errors.

All `(lambda,mu)` irreps are normalized to U(3) `[lambda+mu,mu,0]`.
**All kappa and rho indices in these files are zero-based**, matching the C++
coefficient API. Printed physics notation often uses one-based indices.

## Label schemas

A trailing colon is part of each label.

- `canonical.dat`:
  `A x B -> C : n12A : n22A : n12B : n22B : n12C : n22C : rho :`
  (The coupling string has no spaces.) Values are U(2)-reduced coefficients,
  not full magnetic CGs. Each complete allowed reduced block is recorded.
- `physical.dat`:
  `A x B -> C : LA : LB : LC : kappaA : kappaB : kappaC : rho :`
  Values are SO(3)-reduced coefficients. All allowed L triples and all
  multiplicity indices are included for each documented coupling.
- `conjugation.dat`:
  `(lambda,mu):GT:n12:n22:n11:` stores the Gelfand conjugation sign.
  `(lambda,mu):phase:L:` stores the distinct-irrep phase at M=L;
  the all-M test multiplies by `(-1)^(L-M)`.
  `(n,n):D:L:kappa_target:kappa_source:` stores the full self-conjugate matrix.
- `recoupling.dat`:
  `U6:fundamental:total:left_channel:right_channel:` and likewise `Z6`, `U9`.
  All external factors equal `fundamental` (or the scalar for identity cases).
  For U6 left/right channels are AB/BC; for Z6 they are AB/AC.
  For U9 the left label specifies both AB and CD, the right both AC and BD.
  Every outer multiplicity in these small recoupling examples is one.
- Additional multiplicity-resolved labels in `recoupling.dat`:
  `U6multi:A:B:C:total:AB:rhoAB:rhoAB_C:BC:rhoBC:rhoA_BC:`
  and `Z6multi:A:B:C:total:AB:rhoAB:rhoAB_C:AC:rhoAC:rhoAC_B:`.
  All participating irreps and every zero-based multiplicity index are explicit.
- `U9multi:(1,1):total:AB:CD:rhoAB:rhoCD:rhoAB_CD:AC:BD:rhoAC:rhoBD:rhoAC_BD:`
  records a complete transformation between the `(AB)(CD)` and `(AC)(BD)`
  channel bases for four adjoint irreps. These labels list left then right
  channels, while the structured C++ U9 accessor orders its axes as
  `(rhoAC,rhoBD,rhoAB_CD,rhoAB,rhoCD,rhoAC_BD)`.
- `U9tutorial:A:B:C:D:total:` records the Python tutorial's one-dimensional example;
  AB=CD=AC=BD=(2,0), and all six multiplicity indices are zero.

## Multiplicity recoupling cases

| Coefficients | External irreps | Total | Complete matrix size |
|---|---|---|---|
| U6 and Z6 | (2,0), (0,2), (2,1) | (2,1) | 5 x 5 |
| U6 and Z6 | (1,1), (1,1), (1,1) | (1,1) | 8 x 8 |
| U6 and Z6 | (2,2), (1,1), (1,1) | (2,2) | 10 x 10 |
| U9 | four (1,1) irreps | (0,0) | 8 x 8 |
| U9 | four (1,1) irreps | (1,1) | 32 x 32 |

The first row contains the exact blocks in `python/tutorial.ipynb`, section 10:
AB=(1,1), BC=(1,2) for U6 and AC=(2,2) for Z6. Rather than testing only those
blocks, the suite enumerates every allowed intermediate channel and rho index.
The tutorial's U9 scalar-leg example is included separately as an analytic +1
identity; its multiplicities are all one. The four-adjoint U9 examples add
pair multiplicities, and the total-(1,1) example also has final-coupling
multiplicities. No multiplicity channel is dropped when testing orthogonality.
The suite checks channel counts, coefficient block sizes, and the presence of
nontrivial multiplicities so these examples cannot silently become scalar tests.
All external, intermediate and final irreps in the documented matrices obey
lambda,mu <= 5.

## Conjugation convention

`C|g> = (-1)^(n23-n12-n22-n11) |g.conjugate()>`.
For distinct conjugate irreps:

```
C |(lambda,mu) kappa L M>
  = (-1)^(lambda+mu-L+M) |(mu,lambda) kappa L,-M>.
```

For self-conjugate irreps:

```
C |(n,n) kappa L M>
  = (-1)^(L-M) sum_kappa' D[kappa',kappa] |(n,n) kappa' L,-M>.
```

D is measured in the current library basis, not in a separately diagonalized
basis. The tests do not assume that it is diagonal.

## Analytic checks and recognizable examples

- Coupling a scalar on the left is the identity: the canonical coefficient is
  +1 and the physical coefficient is delta(kappaB,kappaC).
- For R x conjugate(R) -> (0,0), full canonical coefficients are
  `delta(gB, conjugate(gA)) * conjugate_phase(gA) / sqrt(dim(R))`.
  These are checked after restoring the U(2) magnetic CG factor.
- D for (1,1), L=1 is -1; for L=2 it is +1.
- D for (2,2), L=2 is diag(+1,-1).
- The previously measured (4,4), L=2 matrix is numerically consistent with
  `[[4,sqrt(33)],[sqrt(33),-4]]/7`. The test explicitly checks this expression;
  this identification is a regression example, not a general analytic proof.
- With channel ordering [(2,0),(0,1)], the three-fundamental U6 matrix is
  `[[1,sqrt(3)],[sqrt(3),-1]]/2`. Z6 is
  `[[-1,sqrt(3)],[sqrt(3),1]]/2`. The four-fundamental U9 matrix with total
  (0,2) equals the latter matrix in this convention. Scalar recouplings are +1.

## Provenance

Generated separately with `generate_references.cpp`, then checked by
running the test executables against the saved values. Generation uses the
library itself; analytic identities and independent contractions supplement
these regression values rather than making them an external reference.

- Date: 2026-09-29.
- Platform: macOS arm64; GNU C++ (Homebrew GCC) 16.2.0.
- Values are backend-independent to the test tolerances: numerical agreement is
  tolerance-based, not bitwise, across Wigner/solver backends. The current
  libwignernj + krans (dense-QR fallback) stack reproduces every saved value.
- Basis fixing: the library's RREF, Gram-Schmidt and leading-nonzero sign rule;
  `ADD_SU3LIB_PHASE_U3SO3` enabled as in the source.
- Reproduction: run `make config` then `make reference-tool` from the
  u3lib root, then run `./bin/test/generate_references OUTPUT_DIRECTORY`.

No production numerical routines were changed to produce these references.

Initial validation: all five executables passed 104 labeled checks against
5,622 saved values. Altering one reference in an isolated temporary directory
caused the corresponding comparison and executable to fail, as intended.

After adding the tutorial and multiplicity recoupling cases, all five executables
passed 113 labeled checks against 7,089 documented values. Existing reference
values were verified unchanged within tolerance before adding the new entries.
