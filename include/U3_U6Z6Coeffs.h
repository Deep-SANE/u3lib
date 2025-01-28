// protection against multiple header file inclusion
#ifndef U3_U6Z6COEFFS_H
#define U3_U6Z6COEFFS_H

#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

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

    for (int w3 = u3AB.n33; w3 <= u3AB.n13; w3++)
    {
      for (int w1 = u3AB.n13; w1 >= u3AB.n33; w1--)
        if (w3 == u3AB.n33 && w1 == u3AB.n13)
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

    for (int w3 = u3BC.n33; w3 <= u3BC.n13; w3++)
    {
      for (int w1 = u3BC.n13; w1 >= u3BC.n33; w1--)
        if (w3 == u3BC.n33 && w1 == u3BC.n13)
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

    for (int w3 = u3AB.n33; w3 <= u3AB.n13; w3++)
    {
      for (int w1 = u3AB.n13; w1 >= u3AB.n33; w1--)
        if (w3 == u3AB.n33 && w1 == u3AB.n13)
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

    for (int w3 = u3AC.n33; w3 <= u3AC.n13; w3++)
    {
      for (int w1 = u3AC.n13; w1 >= u3AC.n33; w1--)
        if (w3 == u3AC.n33 && w1 == u3AC.n13)
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

/*
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> U6(u3::irrep& u3A,
  u3::irrep& u3B, u3::irrep& u3ABC, u3::irrep& u3C, u3::irrep& u3AB, int
  rhoABmax, int rhoAB_Cmax, u3::irrep& u3BC, int rhoBCmax, int rhoA_BCmax)
  {

    std::vector<u3::gelfand> gArrA = u3::gelfandBasis(u3A);
    std::vector<u3::gelfand> gArrB = u3::gelfandBasis(u3B);
    std::vector<u3::gelfand> gArrC = u3::gelfandBasis(u3C);
    std::vector<u3::gelfand> gArrAB = u3::gelfandBasis(u3AB);
    std::vector<u3::gelfand> gArrBC = u3::gelfandBasis(u3BC);
    std::vector<u3::gelfand> gArrABC = u3::gelfandBasis(u3ABC);

    // coupling <gA;gBC|hwABC> rhoA_BC
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  coeffsA_BC(gArrA.size() * gArrBC.size() * gArrABC.size(), rhoA_BCmax);
    coeffsA_BC.setZero();
    u3::HighestWeightCoeffs(gArrA, gArrBC, gArrABC, coeffsA_BC);

    // coupling <gA;gB|gAB> rhoAB
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsAB(gArrA.size()
  * gArrB.size() * gArrAB.size(), rhoABmax); coeffsAB.setZero();

    for (int w3 = u3AB.n33; w3 <= u3AB.n13; w3++)
    {
      for (int w1 = u3AB.n13; w1 >= u3AB.n33; w1--)
      {
        if (w3 == u3AB.n33 && w1 == u3AB.n13)
        {
          u3::HighestWeightCoeffs(gArrA, gArrB, gArrAB, coeffsAB);
        }
        else
        {
          u3::LowerWeightCoeffs(w3, w1, gArrA, gArrB, gArrAB, coeffsAB);
        }
      }
    }

    // coupling <gAB;gC|hwABC> rhoAB_C
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  coeffsAB_C(gArrAB.size() * gArrC.size() * gArrABC.size(), rhoAB_Cmax);
    coeffsAB_C.setZero();
    u3::HighestWeightCoeffs(gArrAB, gArrC, gArrABC, coeffsAB_C);

    // coupling <gB;gC|hwBC> rhoBC
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsBC(gArrB.size()
  * gArrC.size() * gArrBC.size(), rhoBCmax); coeffsBC.setZero();
    u3::HighestWeightCoeffs(gArrB, gArrC, gArrBC, coeffsBC);


    // U6 coefficients
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matA(rhoA_BCmax,
  rhoA_BCmax), matB(rhoA_BCmax, rhoABmax*rhoAB_Cmax*rhoBCmax); matB.setZero();

    int rA_BC = 0;
    for (int i = 0; i < gArrA.size(); i++)
    {
      if ((gArrA[i].e22() + gArrBC[0].e22()) == gArrABC[0].e22() &&
            (gArrA[i].e11() + gArrBC[0].e11()) == gArrABC[0].e11())
      {
        //std::cout << "Hit" << std::endl;
        for (int j = 0; j < rhoA_BCmax; j++)
          matA(rA_BC,j) = coeffsA_BC( i*gArrBC.size()*gArrABC.size(), j );

        for (int gB = 0; gB < gArrB.size(); gB++)
        {
          for (int gC = 0; gC < gArrC.size(); gC++)
          {
            for (int gAB = 0; gAB < gArrAB.size(); gAB++)
            {
              if ( (gArrA[i].e22() + gArrB[gB].e22()) == gArrAB[gAB].e22() &&
  (gArrA[i].e11() + gArrB[gB].e11()) == gArrAB[gAB].e11()
                  && (gArrAB[gAB].e22() + gArrC[gC].e22()) == gArrABC[0].e22()
  && (gArrAB[gAB].e11() + gArrC[gC].e11()) == gArrABC[0].e11()
                  && (gArrB[gB].e22() + gArrC[gC].e22()) == gArrBC[0].e22() &&
  (gArrB[gB].e11() + gArrC[gC].e11()) == gArrBC[0].e11() )
              {
                for (int rAB = 0; rAB < rhoABmax; rAB++)
                {
                  for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
                  {
                    for (int rBC = 0; rBC < rhoBCmax; rBC++)
                    {
                      matB(rA_BC , (rAB*rhoAB_Cmax + rAB_C)*rhoBCmax + rBC)
                        += coeffsAB((i*gArrB.size() + gB)*gArrAB.size() + gAB,
  rAB)
                        * coeffsAB_C((gAB*gArrC.size() + gC)*gArrABC.size(),
  rAB_C)
                        * coeffsBC((gB*gArrC.size() + gC)*gArrBC.size(), rBC);

                    }
                  }
                }
              }
            }
          }
        }
        rA_BC++;
      }

      if (rA_BC == rhoA_BCmax)
        break;
    }

    //std::cout << "A: " << matA << std::endl;
    //std::cout << "B: " << matB << std::endl;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> U =
  inverse(matA)*matB; //matA.colPivHouseholderQr().solve(matB); return U;
  }

  // input by users: u3b, u3a, u3abc, u3c, u3ab, u3ac
  std::vector<double> Z6(u3::irrep& u3b, u3::irrep& u3a, u3::irrep& u3abc,
                         u3::irrep& u3c, u3::irrep& u3ab,
                         std::vector<int>& rhoab, std::vector<int>& rhoab_c,
                         u3::irrep& u3ac, std::vector<int>& rhoac,
                         std::vector<int>& rhoac_b)
  {
    std::vector<double> Z6coeffs;

    u3::gelfand HW(
        u3abc.n13, u3abc.n23, u3abc.n33, u3abc.n13, u3abc.n23,
        u3abc.n13); // highest weight Gelfand of state of the final irrep
    u3::gelfand HW13(u3ac.n13, u3ac.n23, u3ac.n33, u3ac.n13, u3ac.n23,
                     u3ac.n13);

    // coupling <G13;G2|HW> rhoac_b
    std::vector<u3::gelfand> gArr13ac_b, gArr2ac_b,
        gArrac_b; // store basis states
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsac_b;
    u3::HighestWeightCoeffs(u3ac, u3b, HW, gArr13ac_b, gArr2ac_b, gArrac_b,
                            coeffsac_b);
    int rhoac_b_max = coeffsac_b.cols();

    // indices of <HW13;G2|HW>
    std::vector<int> index_HW13G2;
    for (int i = 0; i < gArr13ac_b.size(); i++)
      if (u3::compare(gArr13ac_b[i], HW13))
        index_HW13G2.push_back(i);

    // coupling <G1;G2|G12> rhoab
    std::vector<u3::gelfand> gArr1ab, gArr2ab, gArr12ab; // store basis states
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsab;
    for (int w3 = u3ab.n33; w3 <= u3ab.n13; w3++)
    {
      for (int w1 = u3ab.n13; w1 >= u3ab.n33; w1--)
      {
        if (w3 == u3ab.n33 && w1 == u3ab.n13)
        {
          u3::gelfand HWab(u3ab.n13, u3ab.n23, u3ab.n33, u3ab.n13, u3ab.n23,
                           u3ab.n13);
          u3::HighestWeightCoeffs(u3a, u3b, HWab, gArr1ab, gArr2ab, gArr12ab,
                                  coeffsab);
        }
        else
        {
          u3::LowerWeightCoeffs(u3a, u3b, u3ab, w3, w1, gArr1ab, gArr2ab,
                                gArr12ab, coeffsab);
        }
      }
    }
    int rhoab_max = coeffsab.cols();

    // coupling <G12;G3|HW> rhoab_c
    std::vector<u3::gelfand> gArr12ab_c, gArr3ab_c,
        gArrab_c; // store basis states
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsab_c;
    u3::HighestWeightCoeffs(u3ab, u3c, HW, gArr12ab_c, gArr3ab_c, gArrab_c,
                            coeffsab_c);
    int rhoab_c_max = coeffsab_c.cols();

    // coupling <G1;G3|HW13> rhoac
    std::vector<u3::gelfand> gArr1ac, gArr3ac, gArr13ac; // store basis states
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsac;
    u3::HighestWeightCoeffs(u3a, u3c, HW13, gArr1ac, gArr3ac, gArr13ac,
                            coeffsac);
    int rhoac_max = coeffsac.cols();

    // Z6 coefficients
    for (int rab = 0; rab < rhoab_max; rab++)
    {
      for (int rab_c = 0; rab_c < rhoab_c_max; rab_c++)
      {
        for (int rac = 0; rac < rhoac_max; rac++)
        {

          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matA(
              rhoac_b_max, rhoac_b_max);
          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matB(
              rhoac_b_max, 1);

          for (int i = 0; i < rhoac_b_max; i++)
          {

            for (int j = 0; j < rhoac_b_max; j++)
              matA(i, j) = coeffsac_b(index_HW13G2[i], j);

            matB(i, 0) = 0.0;
            for (int j = 0; j < coeffsab.rows(); j++)
            {
              if (u3::compare(gArr2ab[j], gArr2ac_b[index_HW13G2[i]]))
              {
                for (int k = 0; k < coeffsab_c.rows(); k++)
                {
                  if (u3::compare(gArr12ab[j], gArr12ab_c[k]))
                  {
                    for (int l = 0; l < coeffsac.rows(); l++)
                    {
                      if (u3::compare(gArr1ac[l], gArr1ab[j]) &&
                          u3::compare(gArr3ac[l], gArr3ab_c[k]))
                      {
                        matB(i, 0) += coeffsab(j, rab) * coeffsab_c(k, rab_c) *
                                      coeffsac(l, rac);
                      }
                    }
                  }
                }
              }
            }
          }

          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Z =
              matA.colPivHouseholderQr().solve(matB);

          for (int rac_b = 0; rac_b < rhoac_b_max; rac_b++)
          {
            rhoab.push_back(rab + 1);
            rhoab_c.push_back(rab_c + 1);
            rhoac.push_back(rac + 1);
            rhoac_b.push_back(rac_b + 1);
            Z6coeffs.push_back(Z(rac_b, 0));
          }
        }
      }
    }
    return Z6coeffs;
  }
*/