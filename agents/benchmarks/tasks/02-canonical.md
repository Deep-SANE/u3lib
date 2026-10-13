# Task 2: Python canonical coefficient batch adapter

Implement `consumer/canonical.py` exposing `canonical_rows(requests)`.
Each request is a tuple (a,b,c); each state is the five integers
(n13,n23,n33,n12,n22), with normalized top row and valid interlacing.
Return a list of independent owning NumPy float64 arrays, one per request,
indexed by outer multiplicity rho. These are U(3)>U(2) reduced coefficients,
not full magnetic Clebsch–Gordan coefficients.

Support input and output irreps with lambda,mu <= 5. Requests satisfy the
subgroup coupling selection rules and have nonzero outer multiplicity. Include
zero coefficient entries and every rho. Returned arrays must retain their
values across subsequent calls and modification of a different returned row.
Repeated calls and an empty batch must work in a fresh process without caller
initialization. The adapter owns the library lifecycle for each call; no other
library client is active. Supply a small runnable example with multiplicity > 1.
