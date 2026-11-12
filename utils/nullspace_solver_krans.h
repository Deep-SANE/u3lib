#ifndef UTILS_NULLSPACE_SOLVER_KRANS_H
#define UTILS_NULLSPACE_SOLVER_KRANS_H

#ifdef USE_KRANS

#include <algorithm>
#include <cstddef>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <krans/Solvers/LOBPCGSolver.hpp>
#include <krans/adapters/eigen/EigenSparseMatrixAdapter.hpp>
#include <krans/adapters/eigen/EigenVectorAdapter.hpp>
#include <krans/adapters/eigen/EigenVectorFactory.hpp>

#include "nullspace_error.h"
#include "nullspace_solver_qr.h" // null_space (QR) + RREF/GramSchmidt/SignConvention

namespace utils
{
  // Right null space of EA via the matrix-free block LOBPCG solver from libkrans
  // (krans). Drop-in sibling of null_space_slepc: same signature
  // and the same QR fall-backs, but needs no PETSc/SLEPc.
  //
  // NOTE on accuracy: libkrans accepts a null vector when ||A v|| <= max(tol_abs,
  // tol_rel * ||A||_F) -- a *relative* criterion -- whereas the post-check below
  // (and the rest of the utils layer) uses an *absolute* threshold of 1e-6. For
  // well-scaled matrices these agree; for large-norm matrices the absolute check
  // may trip and route to the QR fall-back. Revisit the tolerance coupling once
  // benchmarked on real problems.
  inline Eigen::MatrixXd null_space_krans(const Eigen::SparseMatrix<double>& EA,
                                        int nullity)
  {
    const auto n = EA.cols();

    // Mirror null_space_slepc's QR fall-backs for small / block-heavy cases:
    if (EA.size() < 10000) return null_space(EA);
    if (n <= 5 * nullity) return null_space(EA);

    // Matrix-free LOBPCG on N = EA^T EA. libkrans keeps no global state and
    // allocates per-call vectors, so (unlike the SLEPc backend) concurrent
    // calls are safe and no omp critical section is required.
    krans::EigenSparseMatrixAdapter A(EA);
    krans::EigenVectorFactory vf(static_cast<std::size_t>(n));

    krans::solvers::LOBPCGOptions opt;
    opt.nev = static_cast<std::size_t>(nullity);
    opt.nullity_known = nullity; // stop exactly at the target
    opt.block_size =
        std::max<std::size_t>(static_cast<std::size_t>(nullity) + 8, 16);
    opt.tol_abs = 1e-14;
    opt.tol_rel = 1e-14;

    auto basis = krans::solvers::find_nullspace_lobpcg(A, vf, opt);

    // Under-converged -> fall back to the exact QR solver.
    if (static_cast<int>(basis.size()) < nullity) return null_space(EA);

    Eigen::MatrixXd res(n, nullity);
    for (int i = 0; i < nullity; i++)
      res.col(i) = dynamic_cast<krans::EigenVectorAdapter&>(*basis[i]).view();

    // Same canonicalization as null_space / null_space_slepc so the basis
    // follows the existing convention.
#ifdef FIX_MULTIS_AMBIGUITY
    Eigen::MatrixXd res_t = res.transpose();
    ReducedRowEchelon(res_t);
    res = res_t.transpose();

    GramSchmidt(res);
#endif
    // Anchor the column sign on the first genuine coefficient. The LOBPCG
    // path leaves ~1e-14 noise in positions that are exact zeros in the QR
    // reference, so the default 1e-15 threshold would anchor on that noise and
    // flip whole-column signs; genuine coupling coefficients are >> 1e-10.
    SignConvention(res, 1e-12);

    // Inaccurate result -> fall back to QR (graceful; no abort/assert, since
    // this is an opt-in iterative path).
    if (nullspace_error(EA, res) >= 1e-6) return null_space(EA);

    return res;
  }
} // namespace utils

#endif // USE_KRANS

#endif
