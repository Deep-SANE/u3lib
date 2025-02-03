#ifndef U3_U3LIB_U9COEFFS_H
#define U3_U3LIB_U9COEFFS_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "decomposition.h"
#include "U3basis.h"
#include "U3U2U1xU3U2U1.h"

namespace u3
{

  // input by users: u3a, u3b, u3c, u3d, u3ab, u3cd, u3ac, u3bd, u3abcd
  std::vector<double> U9(u3::irrep& u3A, u3::irrep& u3B, u3::irrep& u3C,
                         u3::irrep& u3D, u3::irrep& u3AB, u3::irrep& u3CD,
                         u3::irrep& u3AC, u3::irrep& u3BD, u3::irrep& u3ABCD,
                         int rhoABmax, int rhoCDmax, int rhoAB_CDmax,
                         int rhoACmax, int rhoBDmax, int rhoAC_BDmax)
  {
    std::vector<u3::gelfand> gArrA = u3::gelfandBasis(u3A);
    std::vector<u3::gelfand> gArrB = u3::gelfandBasis(u3B);
    std::vector<u3::gelfand> gArrC = u3::gelfandBasis(u3C);
    std::vector<u3::gelfand> gArrD = u3::gelfandBasis(u3D);
    std::vector<u3::gelfand> gArrAB = u3::gelfandBasis(u3AB);
    std::vector<u3::gelfand> gArrCD = u3::gelfandBasis(u3CD);
    std::vector<u3::gelfand> gArrAC = u3::gelfandBasis(u3AC);
    std::vector<u3::gelfand> gArrBD = u3::gelfandBasis(u3BD);
    std::vector<u3::gelfand> gArrABCD = u3::gelfandBasis(u3ABCD);

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

    // coupling <gC;gD|gCD> rhoCD
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsCD(
        gArrC.size() * gArrD.size() * gArrCD.size(), rhoCDmax);
    coeffsCD.setZero();

    for (int w3 = u3CD.n33(); w3 <= u3CD.n13(); w3++)
    {
      for (int w1 = u3CD.n13(); w1 >= u3CD.n33(); w1--)
        if (w3 == u3CD.n33() && w1 == u3CD.n13())
          u3::HighestWeightCoeffs(gArrC, gArrD, gArrCD, coeffsCD);
        else
          u3::LowerWeightCoeffs(w3, w1, gArrC, gArrD, gArrCD, coeffsCD);
    }

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

    // coupling <gB;gD|gBD> rhoBD
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsBD(
        gArrB.size() * gArrD.size() * gArrBD.size(), rhoBDmax);
    coeffsBD.setZero();

    for (int w3 = u3BD.n33(); w3 <= u3BD.n13(); w3++)
    {
      for (int w1 = u3BD.n13(); w1 >= u3BD.n33(); w1--)
        if (w3 == u3BD.n33() && w1 == u3BD.n13())
          u3::HighestWeightCoeffs(gArrB, gArrD, gArrBD, coeffsBD);
        else
          u3::LowerWeightCoeffs(w3, w1, gArrB, gArrD, gArrBD, coeffsBD);
    }

    // coupling <gAB;gCD|hwABCD> rhoAB_CD
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsAB_CD(
        gArrAB.size() * gArrCD.size() * gArrABCD.size(), rhoAB_CDmax);
    coeffsAB_CD.setZero();
    u3::HighestWeightCoeffs(gArrAB, gArrCD, gArrABCD, coeffsAB_CD);

    // coupling <gAC;gBD|hwABCD> rhoAC_BD
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsAC_BD(
        gArrAC.size() * gArrBD.size() * gArrABCD.size(), rhoAC_BDmax);
    coeffsAC_BD.setZero();
    u3::HighestWeightCoeffs(gArrAC, gArrBD, gArrABCD, coeffsAC_BD);

    std::vector<double> U(rhoACmax * rhoBDmax * rhoAB_CDmax * rhoABmax *
                          rhoCDmax * rhoAC_BDmax);
    for (int i = 0; i < U.size(); i++)
      U[i] = 0.0;

    for (int gA = 0; gA < gArrA.size(); gA++)
    {
      for (int gB = 0; gB < gArrB.size(); gB++)
      {
        for (int gC = 0; gC < gArrC.size(); gC++)
        {
          for (int gD = 0; gD < gArrD.size(); gD++)
          {
            for (int gAB = 0; gAB < gArrAB.size(); gAB++)
            {
              for (int gCD = 0; gCD < gArrCD.size(); gCD++)
              {
                for (int gAC = 0; gAC < gArrAC.size(); gAC++)
                {
                  for (int gBD = 0; gBD < gArrBD.size(); gBD++)
                  {
                    if ((gArrA[gA].e22() + gArrB[gB].e22()) ==
                            gArrAB[gAB].e22() &&
                        (gArrA[gA].e11() + gArrB[gB].e11()) ==
                            gArrAB[gAB].e11() &&
                        (gArrC[gC].e22() + gArrD[gD].e22()) ==
                            gArrCD[gCD].e22() &&
                        (gArrC[gC].e11() + gArrD[gD].e11()) ==
                            gArrCD[gCD].e11() &&
                        (gArrAB[gAB].e22() + gArrCD[gCD].e22()) ==
                            gArrABCD[0].e22() &&
                        (gArrAB[gAB].e11() + gArrCD[gCD].e11()) ==
                            gArrABCD[0].e11() &&
                        (gArrA[gA].e22() + gArrC[gC].e22()) ==
                            gArrAC[gAC].e22() &&
                        (gArrA[gA].e11() + gArrC[gC].e11()) ==
                            gArrAC[gAC].e11() &&
                        (gArrB[gB].e22() + gArrD[gD].e22()) ==
                            gArrBD[gBD].e22() &&
                        (gArrB[gB].e11() + gArrD[gD].e11()) ==
                            gArrBD[gBD].e11() &&
                        (gArrAC[gAC].e22() + gArrBD[gBD].e22()) ==
                            gArrABCD[0].e22() &&
                        (gArrAC[gAC].e11() + gArrBD[gBD].e11()) ==
                            gArrABCD[0].e11())
                    {
                      for (int rAC = 0; rAC < rhoACmax; rAC++)
                      {
                        for (int rBD = 0; rBD < rhoBDmax; rBD++)
                        {
                          for (int rAB_CD = 0; rAB_CD < rhoAB_CDmax; rAB_CD++)
                          {
                            for (int rAB = 0; rAB < rhoABmax; rAB++)
                            {
                              for (int rCD = 0; rCD < rhoCDmax; rCD++)
                              {
                                for (int rAC_BD = 0; rAC_BD < rhoAC_BDmax;
                                     rAC_BD++)
                                {
                                  U[((((rAC * rhoBDmax + rBD) * rhoAB_CDmax +
                                       rAB_CD) *
                                          rhoABmax +
                                      rAB) *
                                         rhoCDmax +
                                     rCD) *
                                        rhoAC_BDmax +
                                    rAC_BD] +=
                                      coeffsAB((gA * gArrB.size() + gB) *
                                                       gArrAB.size() +
                                                   gAB,
                                               rAB) *
                                      coeffsCD((gC * gArrD.size() + gD) *
                                                       gArrCD.size() +
                                                   gCD,
                                               rCD) *
                                      coeffsAB_CD((gAB * gArrCD.size() + gCD) *
                                                      gArrABCD.size(),
                                                  rAB_CD) *
                                      coeffsAC((gA * gArrC.size() + gC) *
                                                       gArrAC.size() +
                                                   gAC,
                                               rAC) *
                                      coeffsBD((gB * gArrD.size() + gD) *
                                                       gArrBD.size() +
                                                   gBD,
                                               rBD) *
                                      coeffsAC_BD((gAC * gArrBD.size() + gBD) *
                                                      gArrABCD.size(),
                                                  rAC_BD);
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
          }
        }
      }
    }

    return U;
  }

} // namespace u3

#endif
