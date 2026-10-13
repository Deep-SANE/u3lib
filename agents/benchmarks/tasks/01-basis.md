# Task 1: C++ basis inventory

Build `consumer/basis_inventory`. Invocation: `basis_inventory lambda mu`.
Support every integer 0 <= lambda,mu <= 5. Print whitespace-delimited lines:
`dimension d`, followed by `branch L multiplicity` for every allowed L in
ascending order, followed by `gelfand n13 n23 n33 n12 n22 n11` for every Gelfand
state in library order. Include every state once; omit absent L values.

This feeds a basis allocation routine: branch multiplicity counts independent
copies of L, not its magnetic degeneracy. Add a local check that the magnetic
state count reconstructed from the branching equals dimension. The grader will
check labels and order, not just the total. Example invocation: `basis_inventory 3 1`.
No coefficient calculations are required.
