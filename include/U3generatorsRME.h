#ifndef U3GENERATORSME_H
#define U3GENERATORSME_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "U3basis.h"
#include "U3SO3_to_U3U2U1.h"

namespace u3
{
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  RME_L(u3::physicalvector& bra, u3::physicalvector& ket,
        std::vector<int>& bra_kappa, std::vector<int>& ket_kappa)
  {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> RME(bra.kappa,
                                                              ket.kappa);

    // running through the multiplicity labels
    for (int i = 1; i <= bra.kappa; i++)
      bra_kappa.push_back(i);

    for (int j = 1; j <= ket.kappa; j++)
      ket_kappa.push_back(j);

    for (int i = 0; i < RME.rows(); i++)
    {
      for (int j = 0; j < RME.cols(); j++)
      {
        if (bra.L != ket.L)
        {
          RME(i, j) = 0.0;
        }
        else
        {
          // Note that my convention for reduced matrix elements follows a
          // standard convention, which involves an extra factor of sqrt(2L+1)
          RME(i, j) = std::sqrt(bra.L * (bra.L + 1.0) * (2.0 * bra.L + 1.0));
        }
      }
    }

    return RME;
  }

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  RME_Q(u3::physicalvector& bra, u3::physicalvector& ket,
        std::vector<int>& bra_kappa, std::vector<int>& ket_kappa)
  {
    // selection rule reads |L'-L| <= 2 <= L'+L
    if (std::abs(bra.L - ket.L) > 2 || (ket.L + bra.L) < 2)
    {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> RME(bra.kappa,
                                                                ket.kappa);

      // running through the multiplicity labels
      for (int i = 1; i <= bra.kappa; i++)
        bra_kappa.push_back(i);

      for (int j = 1; j <= ket.kappa; j++)
        ket_kappa.push_back(j);

      for (int i = 0; i < RME.rows(); i++)
        for (int j = 0; j < RME.cols(); j++)
          RME(i, j) = 0.0;

      return RME;
    }
    else if (bra.L < ket.L)
    {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> RME =
          RME_Q(ket, bra, ket_kappa, bra_kappa);

      int phase = ((bra.L - ket.L) % 2 == 0) ? 1 : -1;

      // running through the multiplicity labels
      for (int i = 0; i < RME.rows(); i++)
        for (int j = 0; j < RME.cols(); j++)
          RME(i, j) *= phase;

      return RME.transpose();
    }
    else
    {
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> RME(bra.kappa,
                                                                ket.kappa);

      // Getting the tensor component of Q
      int m = bra.L - ket.L;

      // running through the multiplicity labels
      for (int i = 1; i <= bra.kappa; i++)
        bra_kappa.push_back(i);

      for (int j = 1; j <= ket.kappa; j++)
        ket_kappa.push_back(j);

      for (int i = 0; i < RME.rows(); i++)
        for (int j = 0; j < RME.cols(); j++)
          RME(i, j) = 0.0;

      // Implementation here are based on Eqs. (51)-(53) in F. Pan et al, Nucl.
      // Phys. A 952, 70-99 (2016) Note that my convention for reduced matrix
      // elements follows a more standard convention, which involves an extra
      // factor of sqrt(2L'+1)
      if (m == 2)
      {
        // Expanding in canonical basis
        std::vector<u3::gelfand> bra_canobasis = u3::canonicalbasis(bra);
        std::vector<u3::gelfand> ket_canobasis = u3::canonicalbasis(ket);

        // Building the projection matrix
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bra_Lproj =
            Lprojection(bra, bra_canobasis);

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ket_Lproj =
            Lprojection(ket, ket_canobasis);

        // and solving its null space
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bra_coeffs =
            phys_to_cano(bra, bra_Lproj);

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ket_coeffs =
            phys_to_cano(ket, ket_Lproj);

        // running through the gelfand states
        for (int k = 0; k < ket_canobasis.size(); k++)
        {
          for (int b = 0; b < bra_canobasis.size(); b++)
          {
            if (bra_canobasis[b].n12 == ket_canobasis[k].n12 &&
                bra_canobasis[b].n22 == ket_canobasis[k].n22 &&
                bra_canobasis[b].n11 == (ket_canobasis[k].n11 + 1))
            {
              // running through the multiplicity labels
              for (int i = 0; i < bra.kappa; i++)
                for (int j = 0; j < ket.kappa; j++)
                  RME(i, j) += std::sqrt(2.0 * bra.L + 1.0) * bra_coeffs(b, i) *
                               ket_coeffs(k, j) * ket_canobasis[k].e12();
            }
          }
        }
      }
      else if (m == 1)
      {
        // Expanding in canonical basis
        std::vector<u3::gelfand> bra_canobasis = u3::canonicalbasis(bra);
        std::vector<u3::gelfand> ket_canobasis = u3::canonicalbasis(ket);

        // Building the projection matrix
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bra_Lproj =
            Lprojection(bra, bra_canobasis);

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ket_Lproj =
            Lprojection(ket, ket_canobasis);

        // and solving its null space
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> bra_coeffs =
            phys_to_cano(bra, bra_Lproj);

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ket_coeffs =
            phys_to_cano(ket, ket_Lproj);

        // running through the gelfand states
        for (int k = 0; k < ket_canobasis.size(); k++)
        {
          for (int b = 0; b < bra_canobasis.size(); b++)
          {
            if (bra_canobasis[b].n12 == (ket_canobasis[k].n12 - 1) &&
                bra_canobasis[b].n22 == ket_canobasis[k].n22 &&
                bra_canobasis[b].n11 == ket_canobasis[k].n11)
            {
              // running through the multiplicity labels
              for (int i = 0; i < bra.kappa; i++)
                for (int j = 0; j < ket.kappa; j++)
                  RME(i, j) += std::sqrt((ket.L + 2.0) * (2.0 * bra.L + 1.0) /
                                         (2.0 * ket.L)) *
                               bra_coeffs(b, i) * ket_coeffs(k, j) *
                               ket_canobasis[k].e32_12();
            }
            else if (bra_canobasis[b].n12 == ket_canobasis[k].n12 &&
                     bra_canobasis[b].n22 == (ket_canobasis[k].n22 - 1) &&
                     bra_canobasis[b].n11 == ket_canobasis[k].n11)
            {
              // running through the multiplicity labels
              for (int i = 0; i < bra.kappa; i++)
                for (int j = 0; j < ket.kappa; j++)
                  RME(i, j) += std::sqrt((ket.L + 2.0) * (2.0 * bra.L + 1.0) /
                                         (2.0 * ket.L)) *
                               bra_coeffs(b, i) * ket_coeffs(k, j) *
                               ket_canobasis[k].e32_22();
            }
            else if (bra_canobasis[b].n12 == ket_canobasis[k].n12 &&
                     bra_canobasis[b].n22 == (ket_canobasis[k].n22 + 1) &&
                     bra_canobasis[b].n11 == (ket_canobasis[k].n11 + 1))
            {
              // running through the multiplicity labels
              for (int i = 0; i < bra.kappa; i++)
                for (int j = 0; j < ket.kappa; j++)
                  RME(i, j) -= std::sqrt((ket.L + 2.0) * (2.0 * bra.L + 1.0) /
                                         (2.0 * ket.L)) *
                               bra_coeffs(b, i) * ket_coeffs(k, j) *
                               ket_canobasis[k].e13_22();
            }
            else if (bra_canobasis[b].n12 == (ket_canobasis[k].n12 + 1) &&
                     bra_canobasis[b].n22 == ket_canobasis[k].n22 &&
                     bra_canobasis[b].n11 == (ket_canobasis[k].n11 + 1))
            {
              // running through the multiplicity labels
              for (int i = 0; i < bra.kappa; i++)
                for (int j = 0; j < ket.kappa; j++)
                  RME(i, j) -= std::sqrt((ket.L + 2.0) * (2 * bra.L + 1.0) /
                                         (2.0 * ket.L)) *
                               bra_coeffs(b, i) * ket_coeffs(k, j) *
                               ket_canobasis[k].e13_12();
            }
          }
        }
      }
      else if (m == 0)
      {
        // Expanding in canonical basis
        std::vector<u3::gelfand> ket_canobasis =
            u3::canonicalbasis(ket); // only need to expand one state because
                                     // bra and ket have the same L

        // Building the projection matrix
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ket_Lproj =
            Lprojection(ket, ket_canobasis);

        // and solving its null space
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ket_coeffs =
            phys_to_cano(ket, ket_Lproj);

        // bra and ket have the same kappa_max
        for (int i = 0; i < bra.kappa; i++)
        {
          for (int j = 0; j < ket.kappa; j++)
          {
            // Eq(53)
            for (int k = 0; k < ket_canobasis.size(); k++)
            {
              // std::cout << ket_canobasis[k] << std::endl;
              RME(i, j) += (ket_canobasis[k].n11 - ket.L) * ket_coeffs(k, i) *
                           ket_coeffs(k, j);
            }
            RME(i, j) *= 6.0;

            if (i == j)
              RME(i, j) +=
                  (3.0 * ket.L - 2.0 * ket.toprow.n13 - 2.0 * ket.toprow.n23);

            RME(i, j) *= std::sqrt((ket.L + 1.0) * (2.0 * ket.L + 1.0) *
                                   (2.0 * ket.L + 3.0) /
                                   (6.0 * ket.L * (2.0 * ket.L - 1.0)));
          }
        }
      }

      return RME;
    }
  }

} // namespace u3

#endif
