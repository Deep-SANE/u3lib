Recommended improvements:

1. utils/inverseSolver.h: I set the sparsity boundary based on my experience with the nullspace solver. I have not actually tested what a good sparsity should be to switch from dense-matrix solver to sparse-matrix solver.

2. utils/nullspaceSolver.h: I set the sparsity boundary based on my experience. However, I remember from our conversation during last summer that sparsity should be much higher for a matrix to be consider as sparse. 

3. Storage of canonical coupling coefficients in program U3U2U1xU3U2U1.h
- Two irreps u3A and u3B coupled to u3C. The gelfand states assosicated with irreps are stored arrays called gArrA, gArrB, gArrC.
- The canonical coupling coefficients < gA; gB | gC > are stord in a matrix of size (gArrA.size() * gArrB.size() * gArrC.size(), rhomax)
- Previously, I store gA, gB, gC explicitly. Now, I only store there indices iA, iB, iC within gArrA, gArrB, gArrC.
That is, each < gA; gB | gC > is mapped to the ((iA*gArrB.size() + iB)*gArrC.size() + iC)-th row of the coefficient matrix. (The column feature is related to outer multiplicity rho).
- This size of the matrix is based on the assumption that every pair of (gA,gB) could couple to a state gC. However, this is not the case, there is a selection rule that limit which (gA,gB) pair can exist for a particular gC. You can see the meaning of this by commenting out the "if" statement on line 74 of src/U3U2U1xU3U2U1.cpp and recompile it. You will see that there are many zeros, which are simply because they are not allowed due to selection rule. 
- To mitigate this, I propose implement a fast lookup table with key to be a tuple of (iA,iB,iC) and storing only a row for each key instead of the entire coefficient matrix.
It could look like this:
	std::unordered_map< std::tuple<int, int, int>, Eigen::Matrix<double, 1, Eigen::Dynamic> > umap;
or like this
	std::unordered_map<std::tuple<int, int, int>, std::vector<double, rhomax> > umap;
- I am not so familiar with fast lookup table and hash function, so I hope Daniel and Tomas could help with this. 