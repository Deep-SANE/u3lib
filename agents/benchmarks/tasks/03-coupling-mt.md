# Task 3: C++ OpenMP batch service for canonical and physical coupling

Build `consumer/coupling_batch`. Read requests from stdin until EOF; each line is
one request, tagged by a leading token:

- `C` — canonical U(3)>U(2) reduced coupling. Followed by 15 integers: for each
  state A B C the five labels `(n13,n23,n33,n12,n22)`, with normalized top row
  and valid interlacing.
- `P` — physical SO(3)-reduced coupling. Followed by 9 integers: `(lambda,mu,L)`
  for A B C.

All labels are `<=5`; every request satisfies the subgroup/angular-momentum
selection rules and has positive multiplicities.

**Evaluate the whole batch in parallel with OpenMP.** Open one `#pragma omp
parallel` region; each worker thread calls `u3::initialize_thread()` on entry and
`u3::finalize_thread()` on exit, and the requests are distributed with
`#pragma omp for`. Collect every result into per-request storage during the
parallel region, and only then — on the main thread, in ascending request order
`q` (zero-based) — serialize the output. Parallel evaluation must not change the
output versus a serial run.

For each request `q`, in request order:

- canonical: print `canon q rhomax`, then for `rho = 0 .. rhomax-1` a line
  `c q rho value`. `value` is the `wigU3U2db().get(...)` coefficient (reduced,
  over the outer multiplicity `rho`).
- physical: print `phys q kappaAmax kappaBmax kappaCmax rhomax`, then for every
  `(kappaA,kappaB,kappaC,rho)` in lexicographic index order a line
  `p q kappaA kappaB kappaC rho value`, from `wigU3SO3db().get_U3SO3(...)`.

Use the library's fixed basis, signs and normalization, including exact
selection-rule zeros. An empty stream produces no output and exits successfully.
The program must build under the configured OpenMP flags and run correctly for
any thread count (`OMP_NUM_THREADS`), with no data races and deterministic
output. Supply a local batch test that mixes `C` and `P` requests, includes a
case with multiplicity `> 1`, and runs under at least two threads.
