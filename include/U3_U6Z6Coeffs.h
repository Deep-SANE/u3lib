#ifndef U3LIB_U3_U6Z6COEFFS_H
#define U3LIB_U3_U6Z6COEFFS_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "decomposition.h"
#include "U3basis.h"
#include "U3U2U1xU3U2U1.h"

namespace u3
{
  // input by users: u3A, u3B, u3ABC, u3C, u3AB, u3BC
  std::vector<double> U6(u3::irrep& u3A, u3::irrep& u3B, u3::irrep& u3ABC,
                         u3::irrep& u3C, u3::irrep& u3AB, int rhoABmax,
                         int rhoAB_Cmax, u3::irrep& u3BC, int rhoBCmax,
                         int rhoA_BCmax)
  {
    std::vector<u3::gelfand> gArrA = u3::gelfandBasis(u3A);
    std::vector<u3::gelfand> gArrB = u3::gelfandBasis(u3B);
    std::vector<u3::gelfand> gArrC = u3::gelfandBasis(u3C);
    std::vector<u3::gelfand> gArrAB = u3::gelfandBasis(u3AB);
    std::vector<u3::gelfand> gArrBC = u3::gelfandBasis(u3BC);
    std::vector<u3::gelfand> gArrABC = u3::gelfandBasis(u3ABC);

    // coupling <gA;gB|gAB> rhoAB
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsAB(
        gArrA.size() * gArrB.size() * gArrAB.size(), rhoABmax);
    coeffsAB.setZero();

    for (int w3 = u3AB.n33(); w3 <= u3AB.n13(); w3++)
    {
      for (int w1 = u3AB.n13(); w1 >= u3AB.n33(); w1--)
        if (w3 == u3AB.n33() && w1 == u3AB.n13())
          u3::HighestWeightCoeffs(gArrA, gArrB, gArrAB, coeffsAB);
        else
          u3::LowerWeightCoeffs(w3, w1, gArrA, gArrB, gArrAB, coeffsAB);
    }

    // coupling <gAB;gC|hwABC> rhoAB_C
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsAB_C(
        gArrAB.size() * gArrC.size() * gArrABC.size(), rhoAB_Cmax);
    coeffsAB_C.setZero();
    u3::HighestWeightCoeffs(gArrAB, gArrC, gArrABC, coeffsAB_C);

    // coupling <gB;gC|gBC> rhoBC
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsBC(
        gArrB.size() * gArrC.size() * gArrBC.size(), rhoBCmax);
    coeffsBC.setZero();

    for (int w3 = u3BC.n33(); w3 <= u3BC.n13(); w3++)
    {
      for (int w1 = u3BC.n13(); w1 >= u3BC.n33(); w1--)
        if (w3 == u3BC.n33() && w1 == u3BC.n13())
          u3::HighestWeightCoeffs(gArrB, gArrC, gArrBC, coeffsBC);
        else
          u3::LowerWeightCoeffs(w3, w1, gArrB, gArrC, gArrBC, coeffsBC);
    }

    // coupling <gA;gBC|hwABC> rhoA_BC
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsA_BC(
        gArrA.size() * gArrBC.size() * gArrABC.size(), rhoA_BCmax);
    coeffsA_BC.setZero();
    u3::HighestWeightCoeffs(gArrA, gArrBC, gArrABC, coeffsA_BC);

    // U6 coefficients
    std::vector<double> U(rhoABmax * rhoAB_Cmax * rhoBCmax * rhoA_BCmax);
    for (int i = 0; i < U.size(); i++)
      U[i] = 0.0;

    for (int gA = 0; gA < gArrA.size(); gA++)
    {
      for (int gB = 0; gB < gArrB.size(); gB++)
      {
        for (int gC = 0; gC < gArrC.size(); gC++)
        {
          for (int gAB = 0; gAB < gArrAB.size(); gAB++)
          {
            for (int gBC = 0; gBC < gArrBC.size(); gBC++)
            {
              if ((gArrA[gA].e22() + gArrB[gB].e22()) == gArrAB[gAB].e22() &&
                  (gArrA[gA].e11() + gArrB[gB].e11()) == gArrAB[gAB].e11() &&
                  (gArrAB[gAB].e22() + gArrC[gC].e22()) == gArrABC[0].e22() &&
                  (gArrAB[gAB].e11() + gArrC[gC].e11()) == gArrABC[0].e11() &&
                  (gArrB[gB].e22() + gArrC[gC].e22()) == gArrBC[gBC].e22() &&
                  (gArrB[gB].e11() + gArrC[gC].e11()) == gArrBC[gBC].e11() &&
                  (gArrA[gA].e22() + gArrBC[gBC].e22()) == gArrABC[0].e22() &&
                  (gArrA[gA].e11() + gArrBC[gBC].e11()) == gArrABC[0].e11())
              {
                for (int rAB = 0; rAB < rhoABmax; rAB++)
                {
                  for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
                  {
                    for (int rBC = 0; rBC < rhoBCmax; rBC++)
                    {
                      for (int rA_BC = 0; rA_BC < rhoA_BCmax; rA_BC++)
                      {
                        U[((rAB * rhoAB_Cmax + rAB_C) * rhoBCmax + rBC) *
                              rhoA_BCmax +
                          rA_BC] +=
                            coeffsAB((gA * gArrB.size() + gB) * gArrAB.size() +
                                         gAB,
                                     rAB) *
                            coeffsAB_C((gAB * gArrC.size() + gC) *
                                           gArrABC.size(),
                                       rAB_C) *
                            coeffsBC((gB * gArrC.size() + gC) * gArrBC.size() +
                                         gBC,
                                     rBC) *
                            coeffsA_BC((gA * gArrBC.size() + gBC) *
                                           gArrABC.size(),
                                       rA_BC);
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

    return U;
  }

  // input by users: u3B, u3B, u3ABC, u3C, u3AB, u3AC
  std::vector<double> Z6(u3::irrep& u3B, u3::irrep& u3A, u3::irrep& u3ABC,
                         u3::irrep& u3C, u3::irrep& u3AB, int rhoABmax,
                         int rhoAB_Cmax, u3::irrep& u3AC, int rhoACmax,
                         int rhoAC_Bmax)
  {
    std::vector<u3::gelfand> gArrA = u3::gelfandBasis(u3A);
    std::vector<u3::gelfand> gArrB = u3::gelfandBasis(u3B);
    std::vector<u3::gelfand> gArrC = u3::gelfandBasis(u3C);
    std::vector<u3::gelfand> gArrAB = u3::gelfandBasis(u3AB);
    std::vector<u3::gelfand> gArrAC = u3::gelfandBasis(u3AC);
    std::vector<u3::gelfand> gArrABC = u3::gelfandBasis(u3ABC);

    // coupling <gA;gB|gAB> rhoAB
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsAB(
        gArrA.size() * gArrB.size() * gArrAB.size(), rhoABmax);
    coeffsAB.setZero();

    for (int w3 = u3AB.n33(); w3 <= u3AB.n13(); w3++)
    {
      for (int w1 = u3AB.n13(); w1 >= u3AB.n33(); w1--)
        if (w3 == u3AB.n33() && w1 == u3AB.n13())
          u3::HighestWeightCoeffs(gArrA, gArrB, gArrAB, coeffsAB);
        else
          u3::LowerWeightCoeffs(w3, w1, gArrA, gArrB, gArrAB, coeffsAB);
    }

    // coupling <gAB;gC|hwABC> rhoAB_C
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsAB_C(
        gArrAB.size() * gArrC.size() * gArrABC.size(), rhoAB_Cmax);
    coeffsAB_C.setZero();
    u3::HighestWeightCoeffs(gArrAB, gArrC, gArrABC, coeffsAB_C);

    // coupling <gA;gC|gAC> rhoAC
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsAC(
        gArrA.size() * gArrC.size() * gArrAC.size(), rhoACmax);
    coeffsAC.setZero();

    for (int w3 = u3AC.n33(); w3 <= u3AC.n13(); w3++)
    {
      for (int w1 = u3AC.n13(); w1 >= u3AC.n33(); w1--)
        if (w3 == u3AC.n33() && w1 == u3AC.n13())
          u3::HighestWeightCoeffs(gArrA, gArrC, gArrAC, coeffsAC);
        else
          u3::LowerWeightCoeffs(w3, w1, gArrA, gArrC, gArrAC, coeffsAC);
    }

    // coupling <gA;gBC|hwABC> rhoAC_B
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsAC_B(
        gArrAC.size() * gArrB.size() * gArrABC.size(), rhoAC_Bmax);
    coeffsAC_B.setZero();
    u3::HighestWeightCoeffs(gArrAC, gArrB, gArrABC, coeffsAC_B);

    // U6 coefficients
    std::vector<double> Z(rhoABmax * rhoAB_Cmax * rhoACmax * rhoAC_Bmax);
    for (int i = 0; i < Z.size(); i++)
      Z[i] = 0.0;

    for (int gA = 0; gA < gArrA.size(); gA++)
    {
      for (int gB = 0; gB < gArrB.size(); gB++)
      {
        for (int gC = 0; gC < gArrC.size(); gC++)
        {
          for (int gAB = 0; gAB < gArrAB.size(); gAB++)
          {
            for (int gAC = 0; gAC < gArrAC.size(); gAC++)
            {
              if ((gArrA[gA].e22() + gArrB[gB].e22()) == gArrAB[gAB].e22() &&
                  (gArrA[gA].e11() + gArrB[gB].e11()) == gArrAB[gAB].e11() &&
                  (gArrAB[gAB].e22() + gArrC[gC].e22()) == gArrABC[0].e22() &&
                  (gArrAB[gAB].e11() + gArrC[gC].e11()) == gArrABC[0].e11() &&
                  (gArrA[gA].e22() + gArrC[gC].e22()) == gArrAC[gAC].e22() &&
                  (gArrA[gA].e11() + gArrC[gC].e11()) == gArrAC[gAC].e11() &&
                  (gArrAC[gAC].e22() + gArrB[gB].e22()) == gArrABC[0].e22() &&
                  (gArrAC[gAC].e11() + gArrB[gB].e11()) == gArrABC[0].e11())
              {
                for (int rAB = 0; rAB < rhoABmax; rAB++)
                {
                  for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
                  {
                    for (int rAC = 0; rAC < rhoACmax; rAC++)
                    {
                      for (int rAC_B = 0; rAC_B < rhoAC_Bmax; rAC_B++)
                      {
                        Z[((rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC) *
                              rhoAC_Bmax +
                          rAC_B] +=
                            coeffsAB((gA * gArrB.size() + gB) * gArrAB.size() +
                                         gAB,
                                     rAB) *
                            coeffsAB_C((gAB * gArrC.size() + gC) *
                                           gArrABC.size(),
                                       rAB_C) *
                            coeffsAC((gA * gArrC.size() + gC) * gArrAC.size() +
                                         gAC,
                                     rAC) *
                            coeffsAC_B((gAC * gArrB.size() + gB) *
                                           gArrABC.size(),
                                       rAC_B);
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

    return Z;
  }

} // namespace u3
#endif
