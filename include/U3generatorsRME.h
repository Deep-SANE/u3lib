#ifndef U3LIB_U3GENERATORSRME_H
#define U3LIB_U3GENERATORSRME_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "U3basis.h"
#include "U3SO3_to_U3U2U1.h"

namespace u3
{
  inline Eigen::MatrixXd RME_L(u3::physicalvector& bra, u3::physicalvector& ket,
        std::vector<int>& bra_kappa, std::vector<int>& ket_kappa)
  {
    Eigen::MatrixXd RME(bra.kappa(), ket.kappa());

    // running through the multiplicity labels
    for (int i = 0; i < bra.kappa(); i++)
      bra_kappa.push_back(i);

    for (int j = 0; j < ket.kappa(); j++)
      ket_kappa.push_back(j);

    for (int i = 0; i < RME.rows(); i++)
    {
      for (int j = 0; j < RME.cols(); j++)
      {
        if (bra.L() != ket.L())
        {
          RME(i, j) = 0.0;
        }
        else
        {
          // Note that my convention for reduced matrix elements follows a
          // standard convention, which involves an extra factor of sqrt(2L+1)
          RME(i, j) = std::sqrt(bra.L() * (bra.L() + 1.0) * (2.0 * bra.L() + 1.0));
        }
      }
    }

    return RME;
  }

  inline Eigen::MatrixXd RME_Q(u3::physicalvector& bra, u3::physicalvector& ket,
        std::vector<int>& bra_kappa, std::vector<int>& ket_kappa)
  {
    Eigen::MatrixXd RME(bra.kappa(), ket.kappa());
    RME.setZero();

    // Getting the tensor component of Q
    int m = bra.L() - ket.L();

    // running through the multiplicity labels
    for (int i = 0; i < bra.kappa(); i++)
      bra_kappa.push_back(i);

    for (int j = 0; j < ket.kappa(); j++)
      ket_kappa.push_back(j);

    // Implementation here are based on Eqs. (51)-(53) in F. Pan et al, Nucl.
    // Phys. A 952, 70-99 (2016). Note that my convention for reduced matrix
    // elements follows a more standard convention, which involves an extra
    // factor of sqrt(2L'+1)
    if (m == 2)
    {
      // Expanding in canonical basis
      std::vector<u3::gelfand> bra_gArr = u3::gelfandBasis(bra.irr());
      auto [bra_indices, bra_coeffs] = u3::phys_to_cano(bra, bra.L(), bra_gArr);

      std::vector<u3::gelfand> ket_gArr = u3::gelfandBasis(ket.irr());
      auto [ket_indices, ket_coeffs] = u3::phys_to_cano(ket, ket.L(), ket_gArr);

      double factor = std::sqrt(2.0 * bra.L() + 1.0);

      // running through the gelfand states
      for (int k = 0; k < ket_indices.size(); k++)
      {
        const auto& ket_g = ket_gArr[ket_indices[k]];
        for (int b = 0; b < bra_indices.size(); b++)
        {
          const auto& bra_g = bra_gArr[bra_indices[b]];
           
          if (u3::E12(ket_g, bra_g))
          {
            // running through the multiplicity labels
            for (int i = 0; i < bra.kappa(); i++)
              for (int j = 0; j < ket.kappa(); j++)
                RME(i, j) += factor * bra_coeffs(b, i) 
                           * ket_coeffs(k, j) * ket_g.e12();
          }
        }
      }
    }
    else if (m == 1)
    {
      // Expanding in canonical basis
      std::vector<u3::gelfand> bra_gArr = u3::gelfandBasis(bra.irr());
      auto [bra_indices, bra_coeffs] = u3::phys_to_cano(bra, bra.L(), bra_gArr);

      std::vector<u3::gelfand> ket_gArr = u3::gelfandBasis(ket.irr());
      auto [ket_indices, ket_coeffs] = u3::phys_to_cano(ket, ket.L(), ket_gArr);

      double factor =
          std::sqrt((ket.L() + 2.0) * (2.0 * bra.L() + 1.0) / (2.0 * ket.L()));

      // running through the gelfand states
      for (int k = 0; k < ket_indices.size(); k++)
      {
        const auto& ket_g = ket_gArr[ket_indices[k]];
        for (int b = 0; b < bra_indices.size(); b++)
        {
          const auto& bra_g = bra_gArr[bra_indices[b]];

          if (u3::E23_12(bra_g, ket_g))
          {
            // running through the multiplicity labels
            for (int i = 0; i < bra.kappa(); i++)
              for (int j = 0; j < ket.kappa(); j++)
                RME(i, j) += factor * bra_coeffs(b, i) 
                           * ket_coeffs(k, j) * ket_g.e32_12();
          }
          else if (u3::E23_22(bra_g, ket_g))
          {
            // running through the multiplicity labels
            for (int i = 0; i < bra.kappa(); i++)
              for (int j = 0; j < ket.kappa(); j++)
                RME(i, j) += factor * bra_coeffs(b, i) 
                           * ket_coeffs(k, j) * ket_g.e32_22();
          }
          else if (u3::E13_12(ket_g, bra_g))
          {
            // running through the multiplicity labels
            for (int i = 0; i < bra.kappa(); i++)
              for (int j = 0; j < ket.kappa(); j++)
                RME(i, j) -= factor * bra_coeffs(b, i) 
                           * ket_coeffs(k, j) * ket_g.e13_12();
          }
          else if (u3::E13_22(ket_g, bra_g))
          {
            // running through the multiplicity labels
            for (int i = 0; i < bra.kappa(); i++)
              for (int j = 0; j < ket.kappa(); j++)
                RME(i, j) -= factor * bra_coeffs(b, i) 
                           * ket_coeffs(k, j) * ket_g.e13_22();
          }
        }
      }
    }
    else if (m == 0)
    {
      // Expanding in canonical basis
      std::vector<u3::gelfand> ket_gArr = u3::gelfandBasis(ket.irr());
      auto [ket_indices, ket_coeffs] = u3::phys_to_cano(ket, ket.L(), ket_gArr);

      double factor = std::sqrt((ket.L() + 1.0) * (2.0 * ket.L() + 1.0) 
          * (2.0 * ket.L() + 3.0) / (6.0 * ket.L() * (2.0 * ket.L() - 1.0)));

      // bra and ket have the same kappa_max
      for (int i = 0; i < bra.kappa(); i++)
      {
        for (int j = 0; j < ket.kappa(); j++)
        {
          // Eq(53)
          for (int k = 0; k < ket_indices.size(); k++)
          {
            // std::cout << ket_canobasis[k] << std::endl;
            RME(i, j) += (ket_gArr[ket_indices[k]].n11() - ket.L()) *
                         ket_coeffs(k, i) * ket_coeffs(k, j);
          }
          RME(i, j) *= 6.0;

          if (i == j)
            RME(i, j) +=
                (3.0 * ket.L() - 2.0 * ket.n13() - 2.0 * ket.n23());

          RME(i, j) *= factor;
        }
      }
    }
    else if (m == -1)
    {
      std::vector<int> new_bra_kappa, new_ket_kappa;

      RME = RME_Q(ket, bra, new_ket_kappa, new_bra_kappa);
      int phase = -1;

      // running through the multiplicity labels
      for (int i = 0; i < RME.rows(); i++)
        for (int j = 0; j < RME.cols(); j++)
          RME(i, j) *= phase;
      RME = RME.transpose();
      bra_kappa = new_bra_kappa;
      ket_kappa = new_ket_kappa;
    }
    else if (m == -2)
    {
      std::vector<int> new_bra_kappa, new_ket_kappa;

      RME = RME_Q(ket, bra, new_ket_kappa, new_bra_kappa);
      RME = RME.transpose();
      bra_kappa = new_bra_kappa;
      ket_kappa = new_ket_kappa;
    }

    return RME;
  }

} // namespace u3

#endif
