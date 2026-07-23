#ifndef U3LIB_U3SO3_to_U3U2U1_H
#define U3LIB_U3SO3_to_U3U2U1_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "U3basis.h"

#include <nullspace_solver_qr.h>
#if defined(NULLSPACE_SOLVER_SLEPC)
#include <nullspace_solver_slepc.h>
#elif defined(NULLSPACE_SOLVER_SNS)
#include <nullspace_solver_sns.h>
#endif

namespace u3
{

  // construction of the canonical basis
  inline std::vector<int> 
  expand_to_gelfand(int M, const std::vector<u3::gelfand>& gArr)
  {
    std::vector<int> indices;

    for (int i = 0; i < gArr.size(); i++)
      if (M == gArr[i].e11() - gArr[i].e22())
        indices.push_back(i);

    return indices;
  }

  // finding the transformation by solving a system of linear equations 
  // resulted from acting L+ = sqrt(2)*(E13 + E32)
  inline auto Lprojection(const std::vector<u3::gelfand>& gArr, 
                          const std::vector<int>& indices)
  {
    // ---- single-pass build (two-pass original kept commented below) ----
    // Evaluate each E13/E23 action once: collect nonzero entries and give each
    // connected target row the next contiguous index lazily.
    const int col_num = indices.size();

    std::vector<Eigen::Triplet<double>> ija;
    int row_num = 0;

    for (int i = 0; i < gArr.size(); i++)
    {
      const auto& g_new = gArr[i];
      bool connected = false;

      for (int j = 0; j < col_num; j++)
      {
        const auto& g = gArr[indices[j]];

        double val = 0.0;
        if      (u3::E13_12(g, g_new)) val = g.e13_12();
        else if (u3::E13_22(g, g_new)) val = g.e13_22();
        else if (u3::E23_12(g_new, g)) val = g.e32_12();
        else if (u3::E23_22(g_new, g)) val = g.e32_22();
        else continue;

        connected = true;
        ija.emplace_back(row_num, j, val);
      }

      if (connected) row_num++;
    }

#ifdef NULLSPACE_SPARSE
    Eigen::SparseMatrix<double> projector(row_num, col_num);
    projector.setFromTriplets(ija.begin(), ija.end());
#else
    Eigen::MatrixXd projector = Eigen::MatrixXd::Zero(row_num, col_num);
    for (const auto& t : ija) projector(t.row(), t.col()) = t.value();
#endif

    return projector;

/* ---- original two-pass version (kept for reference) ----
    std::vector<int> new_indices;

    for (int i = 0; i < gArr.size(); i++)
    {
      const auto& g_new = gArr[i];

      for (int j = 0; j < indices.size(); j++)
      {
        const auto& g = gArr[indices[j]];

        // action of E13
        if (u3::E13_12(g, g_new))
        {
          new_indices.push_back(i);
          break;
        }
        else if (u3::E13_22(g, g_new))
        {
          new_indices.push_back(i);
          break;
        }

        // action of E32
        else if (u3::E23_12(g_new, g))
        {
          new_indices.push_back(i);
          break;
        }
        else if (u3::E23_22(g_new, g))
        {
          new_indices.push_back(i);
          break;
        }
      }
    }

    int row_num = new_indices.size(), col_num = indices.size();

#ifdef NULLSPACE_SPARSE
    std::vector<Eigen::Triplet<double>> ija;
#else
    Eigen::MatrixXd projector(row_num, col_num);
#endif

    for (int i = 0; i < row_num; i++)
    {
      const auto& g_new = gArr[new_indices[i]];

      for (int j = 0; j < col_num; j++)
      {
        const auto& g = gArr[indices[j]];

        // action of E13
        if (u3::E13_12(g, g_new))
#ifdef NULLSPACE_SPARSE
          ija.emplace_back(i, j, g.e13_12());
#else
          projector(i, j) = g.e13_12();
#endif
        else if (u3::E13_22(g, g_new))
#ifdef NULLSPACE_SPARSE
          ija.emplace_back(i, j, g.e13_22());
#else
          projector(i, j) = g.e13_22();
#endif

        // action of E32
        else if (u3::E23_12(g_new, g))

#ifdef NULLSPACE_SPARSE
          ija.emplace_back(i, j, g.e32_12());
#else
          projector(i, j) = g.e32_12();
#endif          
        else if (u3::E23_22(g_new, g))
#ifdef NULLSPACE_SPARSE
          ija.emplace_back(i, j, g.e32_22());
#else
          projector(i, j) = g.e32_22();
#endif 

        else
#ifdef NULLSPACE_SPARSE
          ;
#else
          projector(i, j) = 0.0;
#endif
      }
    }

#ifdef NULLSPACE_SPARSE
    Eigen::SparseMatrix<double> projector(row_num, col_num);
    projector.setFromTriplets(ija.begin(), ija.end());
#endif

    return projector;
*/
  }

  // action of the lowering operator L- = sqrt(2)*(E31 + E23)
  inline auto lowerM(
    const u3::physicalvector& p, const int& M, // M is for the new state
    const std::vector<u3::gelfand>& gArr,
    const std::vector<int>& canoBasis,
    const Eigen::MatrixXd& transformation)
  {
    // normalize coefficient
    double n = std::sqrt(2.0 / ((p.L() - M) * (p.L() + M + 1.0)));

    // store new expansion
    std::vector<int> canoBasis_new = u3::expand_to_gelfand(M, gArr);

    Eigen::MatrixXd transformation_new(canoBasis_new.size(), p.kappa());
    transformation_new.setZero();

    // basis states are obtained by acting with two generators E31 and E23
    for (int i = 0; i < canoBasis_new.size(); i++)
    {
      const auto& g_new = gArr[canoBasis_new[i]];

      for (int j = 0; j < canoBasis.size(); j++)
      {
        const auto& g = gArr[canoBasis[j]];

        if (u3::E13_12(g_new, g))
          transformation_new.row(i) += transformation.row(j) * g.e31_12();
        else if (u3::E13_22(g_new, g))
          transformation_new.row(i) += transformation.row(j) * g.e31_22();
        else if (u3::E23_12(g, g_new))
          transformation_new.row(i) += transformation.row(j) * g.e23_12();
        else if (u3::E23_22(g, g_new))
          transformation_new.row(i) += transformation.row(j) * g.e23_22();
      }

      transformation_new.row(i) *= n;

      double tol = 1e-12;

      for (auto& x : transformation_new.row(i))
        if (std::fabs(x) < tol) x = 0.0;
    }
    return std::make_pair(canoBasis_new, transformation_new);
  }

    // action of the lowering operator L+ = sqrt(2)*(E13 + E32)
  inline auto raiseM(
    const u3::physicalvector& p, const int& M, // M is for the new state
    const std::vector<u3::gelfand>& gArr, 
    const std::vector<int>& canoBasis,
    const Eigen::MatrixXd& transformation)
  {
    // normalize coefficient
    double n = std::sqrt(2.0 / ((p.L() + M) * (p.L() - M + 1.0)));

    // store new expansion
    std::vector<int> canoBasis_new = u3::expand_to_gelfand(M, gArr);

    Eigen::MatrixXd transformation_new(canoBasis_new.size(), p.kappa());
    transformation_new.setZero();

    // basis states are obtained by acting with two generators E31 and E23
    for (int i = 0; i < canoBasis_new.size(); i++)
    {
      const auto& g_new = gArr[canoBasis_new[i]];

      for (int j = 0; j < canoBasis.size(); j++)
      {
        const auto& g = gArr[canoBasis[j]];

        if (u3::E13_12(g, g_new))
          transformation_new.row(i) += transformation.row(j) * g.e13_12();
        else if (u3::E13_22(g, g_new))
          transformation_new.row(i) += transformation.row(j) * g.e13_22();
        else if (u3::E23_12(g_new, g))
          transformation_new.row(i) += transformation.row(j) * g.e32_12();
        else if (u3::E23_22(g_new, g))
          transformation_new.row(i) += transformation.row(j) * g.e32_22();
      }

      transformation_new.row(i) *= n;

      double tol = 1e-12;

      for (auto& x : transformation_new.row(i))
        if (std::fabs(x) < tol) x = 0.0;
    }
    return std::make_pair(canoBasis_new, transformation_new);
  }

  // get the transformation for M <= L
#define ADD_SU3LIB_PHASE_U3SO3
  inline std::pair<std::vector<int>, Eigen::MatrixXd> 
    phys_to_cano(const u3::physicalvector& p, int M, 
                 const std::vector<u3::gelfand>& gArr)
  {
    std::vector<int> canoBasis = u3::expand_to_gelfand(p.L(), gArr);
    Eigen::MatrixXd transformation(canoBasis.size(), p.kappa());

    // Step 1: Generate the HWS M = L
#ifdef ADD_SU3LIB_PHASE_U3SO3
    if (p.lam() > p.mu())
    {
      // conjugate irrep
      u3::irrep irr_conj = p.irr().conjugate();
      u3::physicalvector p_conj(irr_conj, p.L());
      std::vector<u3::gelfand> gArr_conj = u3::gelfandBasis(irr_conj);

      // get the transformation
      auto [canoBasis_conj, transformation_conj] =
        u3::phys_to_cano(p_conj, p.L(), gArr_conj);

      for (size_t i = 0; i < canoBasis_conj.size(); i++)
      {
        const auto& g_conj = gArr_conj[canoBasis_conj[i]];
        double phase = ((g_conj.n22() + g_conj.n11() - p.L()) % 2 == 0) 
                     ? 1.0 : -1.0;
        for (size_t j = 0; j < canoBasis.size(); j++)
        {
          const auto& g = gArr[canoBasis[j]];
          if (g.n12() == g_conj.n13() - g_conj.n22() &&
              g.n22() == g_conj.n13() - g_conj.n12() &&
              g.n11() == g_conj.n13() - g_conj.n12() - g_conj.n22() + g_conj.n11())
          {
            transformation.row(j) = transformation_conj.row(i) * phase;
            break;
          }
        }
      }
    }
    else
    {
      // Building the projection matrix
      auto Lproj = u3::Lprojection(gArr, canoBasis);
      // and solving its null space
#if defined(NULLSPACE_SOLVER_SNS)
      transformation = qcx_utils::null_space_sns(Lproj, p.kappa());
#elif defined(NULLSPACE_SOLVER_SLEPC)
      transformation = qcx_utils::null_space_slepc(Lproj, p.kappa());
#else
      transformation = qcx_utils::null_space(Lproj);
#endif
    }

#else // not defining ADD_SU3LIB_PHASE_U3SO3
    // Building the projection matrix
    auto Lproj = u3::Lprojection(gArr, canoBasis);
    // and solving its null space
#if defined(NULLSPACE_SOLVER_SNS)
    transformation = qcx_utils::null_space_sns(Lproj, p.kappa());
#elif defined(NULLSPACE_SOLVER_SLEPC)
    transformation = qcx_utils::null_space_slepc(Lproj, p.kappa());
#else
    transformation = qcx_utils::null_space(Lproj);
#endif

#endif

    // ensure the nullspace is solved correctly
    assert(p.kappa() == transformation.cols());

    // Step 2: Lowering to abs(M)
    for (int m = p.L() - 1; m >= std::fabs(M); m--)
    {
      auto [canoBasis_new, transformation_new] =
        lowerM(p, m, gArr, canoBasis, transformation);

      transformation = std::move(transformation_new);
      canoBasis = std::move(canoBasis_new);
    }

    // Step 3: if M < 0, generate by conjugation
    if (M < 0)
    {
      std::vector<int> canoBasis_new = u3::expand_to_gelfand(M, gArr);
      Eigen::MatrixXd transformation_new(canoBasis_new.size(), p.kappa());

      for (int j = 0; j < canoBasis_new.size(); j++)
      {
        const auto& g_new = gArr[canoBasis_new[j]];

        for (int k = 0; k < canoBasis.size(); k++)
        {
          const auto& g = gArr[canoBasis[k]];

          if (g.n12() == g_new.n12() && g.n22() == g_new.n22() &&
                g.n12() + g.n22() - g.n11() == g_new.n11())
          {
            double phase = ((p.lam() + p.L() + g.n22()) %
                            2 == 0) ? 1.0 : -1.0;
            transformation_new.row(j) = phase * transformation.row(k);
            break;
          }
        }
      }

      transformation = std::move(transformation_new);
      canoBasis = std::move(canoBasis_new);
    }

    return std::make_pair(canoBasis,transformation);
  }

  inline auto phys_to_cano_allM(const u3::physicalvector& p, 
                                const std::vector<u3::gelfand>& gArr)
  {
    std::vector<std::vector<int>> canoBasis(2 * p.L() + 1);
    std::vector<Eigen::MatrixXd> transformation(2 * p.L() + 1);
    
    // Step 1: Generate the HWS M = L
    canoBasis[0] = u3::expand_to_gelfand(p.L(), gArr);

#ifdef ADD_SU3LIB_PHASE_U3SO3
    if (p.lam() > p.mu())
    {
      transformation[0] = Eigen::MatrixXd(canoBasis[0].size(), p.kappa());

      // conjugate irrep
      u3::irrep irr_conj = p.irr().conjugate();
      u3::physicalvector p_conj(irr_conj, p.L());
      std::vector<u3::gelfand> gArr_conj = u3::gelfandBasis(irr_conj);

      // get the transformation
      auto [canoBasis_conj, transformation_conj] =
        u3::phys_to_cano(p_conj, p.L(), gArr_conj);

      for (size_t i = 0; i < canoBasis_conj.size(); i++)
      {
        const auto& g_conj = gArr_conj[canoBasis_conj[i]];

        double phase = ((g_conj.n22() + g_conj.n11() - p.L()) % 2 == 0) 
                     ? 1.0 : -1.0;
        for (size_t j = 0; j < canoBasis[0].size(); j++)
        {
          const auto& g = gArr[canoBasis[0][j]];
          if (g.n12() == g_conj.n13() - g_conj.n22() &&
              g.n22() == g_conj.n13() - g_conj.n12() &&
              g.n11() == g_conj.n13() - g_conj.n12() - g_conj.n22() + g_conj.n11())
          {
            transformation[0].row(j) = transformation_conj.row(i) * phase;
            break;
          }
        }
      }
    }
    else
    {
      // Building the projection matrix
      auto Lproj = u3::Lprojection(gArr, canoBasis[0]);
      // and solving its null space
#if defined(NULLSPACE_SOLVER_SNS)
      transformation[0] = qcx_utils::null_space_sns(Lproj, p.kappa());
#elif defined(NULLSPACE_SOLVER_SLEPC)
      transformation[0] = qcx_utils::null_space_slepc(Lproj, p.kappa());
#else
      transformation[0] = qcx_utils::null_space(Lproj);
#endif
    }

#else // not defining ADD_SU3LIB_PHASE_U3SO3
    // Building the projection matrix
    auto Lproj = u3::Lprojection(gArr, canoBasis[0]);
    // and solving its null space
#if defined(NULLSPACE_SOLVER_SNS)
    transformation[0] = qcx_utils::null_space_sns(Lproj, p.kappa());
#elif defined(NULLSPACE_SOLVER_SLEPC)
    transformation[0] = qcx_utils::null_space_slepc(Lproj, p.kappa());
#else
    transformation[0] = qcx_utils::null_space(Lproj);
#endif

#endif

    // ensure the nullspace is solved correctly
    assert(p.kappa() == transformation[0].cols());

    // Step 2: Lowering to M = 0 (lowerM computes the new expansion itself)
    for (int m = p.L() - 1; m >= 0; m--)
    {
      auto x = p.L() - m;
      std::tie(canoBasis[x], transformation[x]) = lowerM(
        p, m, gArr, canoBasis[x - 1], transformation[x - 1]);
    }

    // Step 3: generate M < 0 by conjugation
    for (int m = -1; m >= -p.L(); m--)
    {
      auto x_m = p.L() - m, x_p = p.L() + m;
      canoBasis[x_m] = u3::expand_to_gelfand(m, gArr);
      transformation[x_m] = 
          Eigen::MatrixXd(canoBasis[x_m].size(), p.kappa());

      for (int j = 0; j < canoBasis[x_m].size(); j++)
      {
        const auto& g_new = gArr[canoBasis[x_m][j]];

        for (int k = 0; k < canoBasis[x_p].size(); k++)
        {
          const auto& g = gArr[canoBasis[x_p][k]];

          if (g.n12() == g_new.n12() && g.n22() == g_new.n22() &&
              g.n12() + g.n22() - g.n11() == g_new.n11())
          {
            double phase = ((p.lam() + p.L() + g.n22()) %
                              2 == 0) ? 1.0 : -1.0;
            transformation[x_m].row(j) = phase 
                * transformation[x_p].row(k);
            break;
          }
        }
      }
    }

    return std::make_pair(canoBasis,transformation);
  }

} // namespace u3

#endif

/*
  // finding the transformation by solving a system of linear equations 
  // resulted from acting L+ = sqrt(2)*(E13+E32)
  inline auto Lprojection(u3::physicalvector& physvec, 
    const std::vector<u3::gelfand>& gArr, std::vector<int>& indices)
  {

    std::vector<int> new_indices;

    for (int i = 0; i < gArr.size(); i++)
    {
      for (int j = 0; j < indices.size(); j++)
      {
        // action of E13
        if (u3::E13_12(gArr[indices[j]], gArr[i]))
        {
          new_indices.push_back(i);
          break;
        }
        else if (u3::E13_22(gArr[indices[j]], gArr[i]))
        {
          new_indices.push_back(i);
          break;
        }

        // action of E32
        else if (u3::E23_12(gArr[i], gArr[indices[j]]))
        {
          new_indices.push_back(i);
          break;
        }
        else if (u3::E23_22(gArr[i], gArr[indices[j]]))
        {
          new_indices.push_back(i);
          break;
        }
      }
    }

    int row_num = new_indices.size(), col_num = indices.size();

#ifdef NULLSPACE_SPARSE
    std::vector<Eigen::Triplet<double>> ija;
#else
    Eigen::MatrixXd projector(row_num, col_num);
#endif

    for (int i = 0; i < row_num; i++)
    {
      for (int j = 0; j < col_num; j++)
      {
        // action of E13
        if (u3::E13_12(gArr[indices[j]], gArr[new_indices[i]]))
#ifdef NULLSPACE_SPARSE
          ija.emplace_back(i, j, gArr[indices[j]].e13_12());
#else
          projector(i, j) = gArr[indices[j]].e13_12();
#endif
        else if (u3::E13_22(gArr[indices[j]], gArr[new_indices[i]]))
#ifdef NULLSPACE_SPARSE
          ija.emplace_back(i, j, gArr[indices[j]].e13_22());
#else
          projector(i, j) = gArr[indices[j]].e13_22();
#endif

        // action of E32
        else if (u3::E23_12(gArr[new_indices[i]], gArr[indices[j]]))

#ifdef NULLSPACE_SPARSE
          ija.emplace_back(i, j, gArr[indices[j]].e32_12());
#else
          projector(i, j) = gArr[indices[j]].e32_12();
#endif          
        else if (u3::E23_22(gArr[new_indices[i]], gArr[indices[j]]))
#ifdef NULLSPACE_SPARSE
          ija.emplace_back(i, j, gArr[indices[j]].e32_22());
#else
          projector(i, j) = gArr[indices[j]].e32_22();
#endif

        else
#ifdef NULLSPACE_SPARSE
          ;
#else
          projector(i, j) = 0.0;
#endif
      }
    }

#ifdef NULLSPACE_SPARSE
    Eigen::SparseMatrix<double> projector(row_num, col_num);
    projector.setFromTriplets(ija.begin(), ija.end());
#endif

    return projector;
  }

  // action of the lowering operator L- = sqrt(2)*(E31 + E23)
  inline auto lowerM(
    const u3::physicalvector& p, const int& M, // M is for the new state
    const std::vector<u3::gelfand>& gArr, 
    const std::vector<int>& canoBasis,
    const Eigen::MatrixXd& transformation)
  {
    // normalize coefficient
    double n = std::sqrt(2.0 / ((p.L - M) * (p.L + M + 1.0)));

    // store new expansion
    std::vector<int> canoBasis_new = u3::expand_to_gelfand(M, gArr);

    Eigen::MatrixXd transformation_new(canoBasis_new.size(), p.kappa);
    transformation_new.setZero();

    // basis states are obtained by acting with two generators E31 and E23
    for (int i = 0; i < canoBasis_new.size(); i++)
    {
      const auto& g_new = gArr[canoBasis_new[i]];

      for (int j = 0; j < canoBasis.size(); j++)
      {
        const auto& g = gArr[canoBasis[j]];

        if (u3::E13_12(gArr[canoBasis_new[i]], gArr[canoBasis[j]]))
          transformation_new.row(i) +=
                transformation.row(j) * gArr[canoBasis[j]].e31_12();
        else if (u3::E13_22(gArr[canoBasis_new[i]], gArr[canoBasis[j]]))
          transformation_new.row(i) +=
                transformation.row(j) * gArr[canoBasis[j]].e31_22();
        else if (u3::E23_12(gArr[canoBasis[j]], gArr[canoBasis_new[i]]))
          transformation_new.row(i) +=
                transformation.row(j) * gArr[canoBasis[j]].e23_12();
        else if (u3::E23_22(gArr[canoBasis[j]], gArr[canoBasis_new[i]]))
          transformation_new.row(i) +=
                transformation.row(j) * gArr[canoBasis[j]].e23_22();

      transformation_new.row(i) *= n;

#ifdef NULLSPACE_SPARSE
        double tol = 1e-12;
#else
        double tol = 1e-12;
#endif

      for (auto& x : transformation_new.row(i))
        if (std::fabs(x) < tol) x = 0.0;
    }
    return std::make_pair(canoBasis_new, transformation_new);
  }
*/
