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

    for (int w3 = u3AB.n33; w3 <= u3AB.n13; w3++)
    {
      for (int w1 = u3AB.n13; w1 >= u3AB.n33; w1--)
        if (w3 == u3AB.n33 && w1 == u3AB.n13)
          u3::HighestWeightCoeffs(gArrA, gArrB, gArrAB, coeffsAB);
        else
          u3::LowerWeightCoeffs(w3, w1, gArrA, gArrB, gArrAB, coeffsAB);
    }

    // coupling <gC;gD|gCD> rhoCD
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsCD(
        gArrC.size() * gArrD.size() * gArrCD.size(), rhoCDmax);
    coeffsCD.setZero();

    for (int w3 = u3CD.n33; w3 <= u3CD.n13; w3++)
    {
      for (int w1 = u3CD.n13; w1 >= u3CD.n33; w1--)
        if (w3 == u3CD.n33 && w1 == u3CD.n13)
          u3::HighestWeightCoeffs(gArrC, gArrD, gArrCD, coeffsCD);
        else
          u3::LowerWeightCoeffs(w3, w1, gArrC, gArrD, gArrCD, coeffsCD);
    }

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

    // coupling <gB;gD|gBD> rhoBD
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsBD(
        gArrB.size() * gArrD.size() * gArrBD.size(), rhoBDmax);
    coeffsBD.setZero();

    for (int w3 = u3BD.n33; w3 <= u3BD.n13; w3++)
    {
      for (int w1 = u3BD.n13; w1 >= u3BD.n33; w1--)
        if (w3 == u3BD.n33 && w1 == u3BD.n13)
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

/*
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  U9(u3::irrep& u3a, u3::irrep& u3b, u3::irrep& u3c, u3::irrep& u3d,
     u3::irrep& u3ab, u3::irrep& u3cd, u3::irrep& u3ac, u3::irrep& u3bd,
     u3::irrep& u3abcd, std::vector<int>& rhoab, int rhoab_max,
     std::vector<int>& rhocd, int rhocd_max, std::vector<int>& rhoab_cd,
     int rhoab_cd_max, std::vector<int>& rhoac, int rhoac_max,
     std::vector<int>& rhobd, int rhobd_max, std::vector<int>& rhoac_bd,
     int rhoac_bd_max)
  {
    // a matrix to store U9-coefficients
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> U9coeffs(
        rhoac_max * rhobd_max * rhoab_cd_max,
        rhoab_max * rhocd_max * rhoac_bd_max);

    // couple the I-III and II irreps
    std::vector<int> rhomaxsac_b;
    std::vector<u3::irrep> u3ac_b = u3::decomposition(u3ac, u3b, rhomaxsac_b);

    // couple the I-II and III irreps
    std::vector<int> rhomaxsab_c;
    std::vector<u3::irrep> u3ab_c = u3::decomposition(u3ab, u3c, rhomaxsab_c);

    // an array to store the couple I-II-III irreps
    std::vector<int> rhomaxsabc;
    std::vector<u3::irrep> u3abc;

    for (int i = 0; i < u3ab_c.size(); i++)
    {
      for (int j = 0; j < u3ac_b.size(); j++)
      {
        if (u3ab_c[i].n13 == u3ac_b[j].n13 && u3ab_c[i].n23 == u3ac_b[j].n23 &&
            u3ab_c[i].n33 == u3ac_b[j].n33)
        {
          u3abc.push_back(u3ab_c[i]);
          rhomaxsabc.push_back(rhomaxsab_c[i]);
          break;
        }
      }
    }

    std::vector<std::vector<int>> rhoab_arr, rhocd_arr, rhoac_arr, rhobd_arr,
        rhoab_cd_arr, rhoac_bd_arr, rhoac_b_uarr, rhoac_b_zarr, rhoab_c_zarr,
        rhoab_c_uarr, rhoabc_d_arr1, rhoabc_d_arr2;
    std::vector<std::vector<double>> U6coeffs_arr1, Z6coeffs_arr, U6coeffs_arr2;

    for (int i = 0; i < u3abc.size(); i++)
    {
      std::vector<int> rab, rcd, rac, rbd, rab_cd, rac_bd, rac_b_u, rac_b_z,
          rab_c_z, rab_c_u, rabc_d_1, rabc_d_2;
      std::vector<double> U6coeffs_1, Z6coeffs, U6coeffs_2;

      U6coeffs_1 = u3::U6(u3ac, u3b, u3abcd, u3d, u3abc[i], rac_b_u, rabc_d_1,
                          u3bd, rbd, rac_bd);
      Z6coeffs = u3::Z6(u3b, u3a, u3abc[i], u3c, u3ab, rab, rab_c_z, u3ac, rac,
                        rac_b_z);
      U6coeffs_2 = u3::U6(u3ab, u3c, u3abcd, u3d, u3abc[i], rab_c_u, rabc_d_2,
                          u3cd, rcd, rab_cd);

      rhoab_arr.push_back(rab);
      rhocd_arr.push_back(rcd);
      rhoac_arr.push_back(rac);
      rhobd_arr.push_back(rbd);
      rhoab_cd_arr.push_back(rab_cd);
      rhoac_bd_arr.push_back(rac_bd);
      rhoac_b_uarr.push_back(rac_b_u);
      rhoac_b_zarr.push_back(rac_b_z);
      rhoab_c_zarr.push_back(rab_c_z);
      rhoab_c_uarr.push_back(rab_c_u);
      rhoabc_d_arr1.push_back(rabc_d_1);
      rhoabc_d_arr2.push_back(rabc_d_2);

      U6coeffs_arr1.push_back(U6coeffs_1);
      U6coeffs_arr2.push_back(U6coeffs_2);
      Z6coeffs_arr.push_back(Z6coeffs);
    }

    int row = 0, col = 0;
    for (int rac = 0; rac < rhoac_max; rac++)
    {
      for (int rbd = 0; rbd < rhobd_max; rbd++)
      {
        for (int rab_cd = 0; rab_cd < rhoab_cd_max; rab_cd++)
        {
          rhoac.push_back(rac + 1);
          rhobd.push_back(rbd + 1);
          rhoab_cd.push_back(rab_cd + 1);
        }
      }
    }

    for (int rab = 0; rab < rhoab_max; rab++)
    {
      for (int rcd = 0; rcd < rhocd_max; rcd++)
      {
        for (int rac_bd = 0; rac_bd < rhoac_bd_max; rac_bd++)
        {
          rhoab.push_back(rab + 1);
          rhocd.push_back(rcd + 1);
          rhoac_bd.push_back(rac_bd + 1);
        }
      }
    }

    for (int m = 0; m < rhoab_cd.size(); m++)
    { // m+1 is the multiplicity label for rhoac, rhobd, rhoab_cd
      for (int n = 0; n < rhoac_bd.size(); n++)
      { // n+1 is the multiplicity label for rhoab, rhocd, rhoac_bd
        U9coeffs(m, n) = 0.0;

        for (int i = 0; i < u3abc.size(); i++)
        {
          for (int j = 0; j < U6coeffs_arr1[i].size(); j++)
          {
            for (int k = 0; k < Z6coeffs_arr[i].size(); k++)
            {
              for (int l = 0; l < U6coeffs_arr2[i].size(); l++)
              {
                if (rhoab_arr[i][k] == rhoab[n] &&
                    rhocd_arr[i][l] == rhocd[n] &&
                    rhoac_bd_arr[i][j] == rhoac_bd[n] &&
                    rhoac_arr[i][k] == rhoac[m] &&
                    rhobd_arr[i][j] == rhobd[m] &&
                    rhoab_cd_arr[i][l] == rhoab_cd[m] &&
                    rhoac_b_uarr[i][j] == rhoac_b_zarr[i][k] &&
                    rhoabc_d_arr1[i][j] == rhoabc_d_arr2[i][l] &&
                    rhoab_c_zarr[i][k] == rhoab_c_uarr[i][l])
                {
                  U9coeffs(m, n) += U6coeffs_arr1[i][j] * Z6coeffs_arr[i][k] *
                                    U6coeffs_arr2[i][l];
                }
              }
            }
          }
        }
      }
    }

    return U9coeffs;
  }
*/
