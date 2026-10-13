#ifndef UTILS_NULLSPACE_SOLVER_SNS_H
#define UTILS_NULLSPACE_SOLVER_SNS_H

#ifdef USE_SNS

#include <algorithm>
#include <cstddef>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <sns/Solvers/LOBPCGSolver.hpp>
#include <sns/adapters/eigen/EigenSparseMatrixAdapter.hpp>
#include <sns/adapters/eigen/EigenVectorAdapter.hpp>
#include <sns/adapters/eigen/EigenVectorFactory.hpp>

#include "nullspace_error.h"
#include "nullspace_solver_qr.h" // null_space (QR) + RREF/GramSchmidt/SignConvention

namespace utils
{
  // Right null space of EA via the matrix-free block LOBPCG solver from libsns
  // (dp-sparse-nullspace). Drop-in sibling of null_space_slepc: same signature
  // and the same QR fall-backs, but needs no PETSc/SLEPc.
  //
  // NOTE on accuracy: libsns accepts a null vector when ||A v|| <= max(tol_abs,
  // tol_rel * ||A||_F) -- a *relative* criterion -- whereas the post-check below
  // (and the rest of the utils layer) uses an *absolute* threshold of 1e-6. For
  // well-scaled matrices these agree; for large-norm matrices the absolute check
  // may trip and route to the QR fall-back. Revisit the tolerance coupling once
  // benchmarked on real problems.
  inline Eigen::MatrixXd null_space_sns(const Eigen::SparseMatrix<double>& EA,
                                        int nullity)
  {
    const auto n = EA.cols();

    // Mirror null_space_slepc's QR fall-backs for small / block-heavy cases:
    if (EA.size() < 10000) return null_space(EA);
    if (n <= 5 * nullity) return null_space(EA);

    // Matrix-free LOBPCG on N = EA^T EA. libsns keeps no global state and
    // allocates per-call vectors, so (unlike the SLEPc backend) concurrent
    // calls are safe and no omp critical section is required.
    sns::EigenSparseMatrixAdapter A(EA);
    sns::EigenVectorFactory vf(static_cast<std::size_t>(n));

    sns::solvers::LOBPCGOptions opt;
    opt.nev = static_cast<std::size_t>(nullity);
    opt.nullity_known = nullity; // stop exactly at the target
    opt.block_size =
        std::max<std::size_t>(static_cast<std::size_t>(nullity) + 8, 16);

    auto basis = sns::solvers::find_nullspace_lobpcg(A, vf, opt);

    // Under-converged -> fall back to the exact QR solver.
    if (static_cast<int>(basis.size()) < nullity) return null_space(EA);

    Eigen::MatrixXd res(n, nullity);
    for (int i = 0; i < nullity; i++)
      res.col(i) = dynamic_cast<sns::EigenVectorAdapter&>(*basis[i]).view();

    // Same canonicalization as null_space / null_space_slepc so the basis
    // follows the existing convention.
#ifdef FIX_MULTIS_AMBIGUITY
    Eigen::MatrixXd res_t = res.transpose();
    ReducedRowEchelon(res_t);
    res = res_t.transpose();

    GramSchmidt(res);
#endif
    SignConvention(res);

    // Inaccurate result -> fall back to QR (graceful; no abort/assert, since
    // this is an opt-in iterative path).
    if (nullspace_error(EA, res) >= 1e-6) return null_space(EA);

    return res;
  }
} // namespace utils

#endif // USE_SNS

#endif
