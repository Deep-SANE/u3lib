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

  // input by users: u3a, u3b, u3abc, u3c, u3ab, u3bc
  std::vector<double> U6(u3::irrep& u3a, u3::irrep& u3b, u3::irrep& u3abc,
                         u3::irrep& u3c, u3::irrep& u3ab,
                         std::vector<int>& rhoab, std::vector<int>& rhoab_c,
                         u3::irrep& u3bc, std::vector<int>& rhobc,
                         std::vector<int>& rhoa_bc)
  {
    std::vector<double> U6coeffs;

    u3::gelfand HW(
        u3abc.n13, u3abc.n23, u3abc.n33, u3abc.n13, u3abc.n23,
        u3abc.n13); // highest weight Gelfand of state of the final irrep
    u3::gelfand HW23(u3bc.n13, u3bc.n23, u3bc.n33, u3bc.n13, u3bc.n23,
                     u3bc.n13);

    // coupling <G1;G23|HW> rhoa_bc
    std::vector<u3::gelfand> gArr1a_bc, gArr23a_bc,
        gArra_bc; // store basis states
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsa_bc;
    u3::HighestWeightCoeffs(u3a, u3bc, HW, gArr1a_bc, gArr23a_bc, gArra_bc,
                            coeffsa_bc);
    int rhoa_bc_max = coeffsa_bc.cols();

    // indices of <G1;HW23|HW>
    std::vector<int> index_G1HW23;
    for (int i = 0; i < gArr1a_bc.size(); i++)
      if (u3::compare(gArr23a_bc[i], HW23))
        index_G1HW23.push_back(i);

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

    // coupling <G2;G3|HW23> rhobc
    std::vector<u3::gelfand> gArr2bc, gArr3bc, gArr23bc; // store basis states
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffsbc;
    u3::HighestWeightCoeffs(u3b, u3c, HW23, gArr2bc, gArr3bc, gArr23bc,
                            coeffsbc);
    int rhobc_max = coeffsbc.cols();

    // U6 coefficients
    for (int rab = 0; rab < rhoab_max; rab++)
    {
      for (int rab_c = 0; rab_c < rhoab_c_max; rab_c++)
      {
        for (int rbc = 0; rbc < rhobc_max; rbc++)
        {

          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matA(
              rhoa_bc_max, rhoa_bc_max);
          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matB(
              rhoa_bc_max, 1);

          for (int i = 0; i < rhoa_bc_max; i++)
          {

            for (int j = 0; j < rhoa_bc_max; j++)
              matA(i, j) = coeffsa_bc(index_G1HW23[i], j);

            matB(i, 0) = 0.0;
            for (int j = 0; j < coeffsab.rows(); j++)
            {
              if (u3::compare(gArr1ab[j], gArr1a_bc[index_G1HW23[i]]))
              {
                for (int k = 0; k < coeffsab_c.rows(); k++)
                {
                  if (u3::compare(gArr12ab[j], gArr12ab_c[k]))
                  {
                    for (int l = 0; l < coeffsbc.rows(); l++)
                    {
                      if (u3::compare(gArr2bc[l], gArr2ab[j]) &&
                          u3::compare(gArr3bc[l], gArr3ab_c[k]))
                      {
                        matB(i, 0) += coeffsab(j, rab) * coeffsab_c(k, rab_c) *
                                      coeffsbc(l, rbc);
                      }
                    }
                  }
                }
              }
            }
          }

          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> U =
              matA.colPivHouseholderQr().solve(matB);

          for (int ra_bc = 0; ra_bc < rhoa_bc_max; ra_bc++)
          {
            rhoab.push_back(rab + 1);
            rhoab_c.push_back(rab_c + 1);
            rhobc.push_back(rbc + 1);
            rhoa_bc.push_back(ra_bc + 1);
            U6coeffs.push_back(U(ra_bc, 0));
          }
        }
      }
    }
    return U6coeffs;
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

} // namespace u3
#endif