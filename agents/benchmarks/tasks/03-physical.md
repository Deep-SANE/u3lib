# Task 3: Python physical coefficient service adapter

Implement `consumer/physical.py` exposing `physical_block(a,b,c)`.
Each argument is (lambda,mu,L). Return an owning NumPy float64 array with axes
(kappaA,kappaB,kappaC,rho) in that order, containing SO(3)-reduced coefficients.
Use the library basis and signs, including exact selection-rule zeros.

Labels must be Python integers (booleans are invalid), 0 <= lambda,mu <= 5,
and 0 <= L <= lambda+mu. Wrong length/type/range raises ValueError before any
unsafe library access. Valid-range inputs with an absent L, failed angular
momentum triangle, or forbidden irrep coupling also raise ValueError. Do not
return an empty array for these cases. Allowed couplings with zero numerical
coefficients remain valid. The process must survive invalid requests and then
successfully answer a valid one.

Repeated calls require no caller initialization; the adapter owns its library
lifecycle and no other client is active. Arrays must remain valid after later
calls. Include local tests covering a nontrivial kappa/rho block and recovery
after invalid input. Example allowed request: (2,2,2),(1,1,2),(2,2,2).
