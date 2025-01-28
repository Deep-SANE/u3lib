#ifndef U3LIB_U3SO3xU3SO3_H
#define U3LIB_U3SO3xU3SO3_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "U3basis.h"
#include "U3SO3_to_U3U2U1.h"
#include "U3U2U1xU3U2U1.h"

extern "C"
{
#include <wigxjpf.h>
}

namespace u3
{

  // usual SO(3) Clebsh-Gordan coefficients <L1 M1; L2 M2 | L3 M3>
  double CGC_SO3(int l1, int l2, int l3, int m1, int m2, int m3)
  {
    int phase = ((-l1 + l2 - m3) % 2 == 0) ? 1 : -1;

    wig_table_init(2 * 100, 9);
    wig_temp_init(2 * 100);

    double val3j = wig3jj(2 * l1, 2 * l2, 2 * l3, 2 * m1, 2 * m2, 2 * (-m3));

    wig_temp_free();
    wig_table_free();

    return phase * std::sqrt(2 * l3 + 1) * val3j;
  }

  // generate all Wigner U3->SO(3) coefficients
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  WigU3SO3(u3::irrep& u3A, u3::irrep& u3B, u3::irrep& u3C, int& rhomax,
           int& dimA, int& dimB, int& dimC,
           std::vector<u3::physicalvector>& pArrA, std::vector<int>& iA,
           std::vector<u3::physicalvector>& pArrB, std::vector<int>& iB,
           std::vector<u3::physicalvector>& pArrC, std::vector<int>& iC)
  {

    // find all the canonical coupling coefficients
    std::vector<u3::gelfand> gArrA;
    std::vector<u3::gelfand> gArrB;
    std::vector<u3::gelfand> gArrC;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> canoCoupleCoeffs =
        u3::WigU3U2U1(u3A, u3B, u3C, rhomax, gArrA, gArrB, gArrC);
    // note that u3C and all gelfand states are already normalized here

    // find physical states for the second irrep
    pArrB = u3::physicalBasis(u3B);
    dimB = 0;
    // store transformation, each column corresponds to one kappa
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
        allCanoCoeffsB(pArrB.size());
    // store expansion
    std::vector<std::vector<int>> allCanoBasisB(pArrB.size());
    // get the transformation
    for (int i = 0; i < pArrB.size(); i++)
    {
      iB.push_back(dimB);
      allCanoCoeffsB[i] = u3::phys_to_cano(pArrB[i], allCanoBasisB[i], gArrB);
      dimB += pArrB[i].kappa;
    }

    // find expansion for the third irrep
    pArrC = u3::physicalBasis(u3C);
    dimC = 0;
    // store transformation, each column corresponds to one kappa
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
        allCanoCoeffsC(pArrC.size());
    // store expansion
    std::vector<std::vector<int>> allCanoBasisC(pArrC.size());
    for (int i = 0; i < pArrC.size(); i++)
    {
      iC.push_back(dimC);
      allCanoCoeffsC[i] = u3::phys_to_cano(pArrC[i], allCanoBasisC[i], gArrC);
      dimC += pArrC[i].kappa;
    }

    // find expansion for the first irrep
    pArrA = u3::physicalBasis(u3A);
    dimA = 0;
    for (int i = 0; i < pArrA.size(); i++)
    {
      iA.push_back(dimA);
      dimA += pArrA[i].kappa;
    }

    // compute Wigner coefficients: < (lamda1,mu1),kappa1,L1 ;
    // (lamda2,mu2),kappa2,L2 | (lamda3,mu3),kappa3,L3 >
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffs(
        dimA * dimB * dimC, rhomax);
    coeffs.setZero();

    for (int i = 0; i < pArrA.size(); i++)
    {

      // std::cout << pArrA[i] << std::endl;

      // store transformation, each column corresponds to one kappa
      std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
          allCanoCoeffsA(2 * pArrA[i].L + 1);
      // store expansion
      std::vector<std::vector<int>> allCanoBasisA(2 * pArrA[i].L + 1);

      allCanoCoeffsA[0] = u3::phys_to_cano(pArrA[i], allCanoBasisA[0], gArrA);

      // now going lower from the HW state, remember that gamma = L - M =
      // 0,1,...,2L
      for (int gamma = 1; gamma <= 2 * pArrA[i].L; gamma++)
      {

        // store new expansion
        allCanoBasisA[gamma] = u3::expand_to_gelfand(pArrA[i].L - gamma, gArrA);

        // normalize coefficient
        double n = std::sqrt(2.0 / (gamma * (2 * pArrA[i].L - gamma + 1.0)));

        allCanoCoeffsA[gamma] =
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(
                allCanoBasisA[gamma].size(), pArrA[i].kappa);
        allCanoCoeffsA[gamma].setZero();

        // basis states are obtained by acting with two generators E31 and E23
        for (int m = 0; m < allCanoBasisA[gamma].size(); m++)
        {
          for (int j = 0; j < allCanoBasisA[gamma - 1].size(); j++)
          {
            if (u3::E13_12(gArrA[allCanoBasisA[gamma][m]],
                           gArrA[allCanoBasisA[gamma - 1][j]]))
            {
              for (int k = 0; k < pArrA[i].kappa; k++)
                allCanoCoeffsA[gamma](m, k) +=
                    n * allCanoCoeffsA[gamma - 1](j, k) *
                    gArrA[allCanoBasisA[gamma - 1][j]].e31_12();
            }
            else if (u3::E13_22(gArrA[allCanoBasisA[gamma][m]],
                                gArrA[allCanoBasisA[gamma - 1][j]]))
            {
              for (int k = 0; k < pArrA[i].kappa; k++)
                allCanoCoeffsA[gamma](m, k) +=
                    n * allCanoCoeffsA[gamma - 1](j, k) *
                    gArrA[allCanoBasisA[gamma - 1][j]].e31_22();
            }
            else if (u3::E23_12(gArrA[allCanoBasisA[gamma - 1][j]],
                                gArrA[allCanoBasisA[gamma][m]]))
            {
              for (int k = 0; k < pArrA[i].kappa; k++)
                allCanoCoeffsA[gamma](m, k) +=
                    n * allCanoCoeffsA[gamma - 1](j, k) *
                    gArrA[allCanoBasisA[gamma - 1][j]].e23_12();
            }
            else if (u3::E23_22(gArrA[allCanoBasisA[gamma - 1][j]],
                                gArrA[allCanoBasisA[gamma][m]]))
            {
              for (int k = 0; k < pArrA[i].kappa; k++)
                allCanoCoeffsA[gamma](m, k) +=
                    n * allCanoCoeffsA[gamma - 1](j, k) *
                    gArrA[allCanoBasisA[gamma - 1][j]].e23_22();
            }
          }
        }

        // std::cout << gamma << " ... " << allCanoCoeffsA[gamma] << std::endl;
      }

      // last step is to compute Wigner coefficients
      // < (lamda1,mu1)kappa1 L1; (lamda2,mu2)kappa2 L2 | (lamda3,mu3)kappa3 L3
      // >
      for (int m = 0; m < pArrB.size(); m++)
      {
        for (int n = 0; n < pArrC.size(); n++)
        {
          if (abs(pArrA[i].L - pArrB[m].L) <= pArrC[n].L &&
              (pArrA[i].L + pArrB[m].L) >= pArrC[n].L)
          {
            double so3cgc =
                CGC_SO3(pArrA[i].L, pArrB[m].L, pArrC[n].L,
                        pArrC[n].L - pArrB[m].L, pArrB[m].L, pArrC[n].L);

            int gamma = pArrA[i].L + pArrB[m].L - pArrC[n].L;

            for (int c1 = 0; c1 < allCanoBasisA[gamma].size(); c1++)
            {
              for (int c2 = 0; c2 < allCanoBasisB[m].size(); c2++)
              {
                for (int c3 = 0; c3 < allCanoBasisC[n].size(); c3++)
                {
                  for (int k1 = 0; k1 < pArrA[i].kappa; k1++)
                  {
                    for (int k2 = 0; k2 < pArrB[m].kappa; k2++)
                    {
                      for (int k3 = 0; k3 < pArrC[n].kappa; k3++)
                      {
                        double x = allCanoCoeffsA[gamma](c1, k1) *
                                   allCanoCoeffsB[m](c2, k2) *
                                   allCanoCoeffsC[n](c3, k3) / so3cgc;
                        for (int r = 0; r < rhomax; r++)
                        {
                          // std::cout << ((i+k1)*dimB + m+k2 )*dimC + n+k3 <<
                          // std::endl;
                          coeffs(((iA[i] + k1) * dimB + iB[m] + k2) * dimC +
                                     iC[n] + k3,
                                 r) +=
                              x * canoCoupleCoeffs(
                                      (allCanoBasisA[gamma][c1] * gArrB.size() +
                                       allCanoBasisB[m][c2]) *
                                              gArrC.size() +
                                          allCanoBasisC[n][c3],
                                      r);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // std::cout << coeffs << std::endl;

    return coeffs;
  }

} // namespace u3

#endif
