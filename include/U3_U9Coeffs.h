#ifndef U3LIB_U3_U9COEFFS_H
#define U3LIB_U3_U9COEFFS_H 

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "U3decomposition.h"
#include "U3_U6Z6Coeffs.h"
#include "U3basis.h"

#include <so3Coeffs/SO3_CGC.h>
#include <so3Coeffs/SO3_Racah.h>

namespace u3
{

  // input by users: u3a, u3b, u3c, u3d, u3ab, u3cd, u3ac, u3bd, u3abcd
  // all normalized
  inline std::vector<double> 
    U9(const u3::irrep& u3A, const std::vector<u3::gelfand>& gArrA,
       const u3::irrep& u3B, const std::vector<u3::gelfand>& gArrB,
       const u3::irrep& u3C, const std::vector<u3::gelfand>& gArrC,
       const u3::irrep& u3D, const std::vector<u3::gelfand>& gArrD,
       const u3::irrep& u3AB, const std::vector<u3::gelfand>& gArrAB,
       const u3::irrep& u3CD, const std::vector<u3::gelfand>& gArrCD,
       const u3::irrep& u3AC, const std::vector<u3::gelfand>& gArrAC,
       const u3::irrep& u3BD, const std::vector<u3::gelfand>& gArrBD,
       const u3::irrep& u3ABCD, const std::vector<u3::gelfand>& gArrABCD,
       int rhoABmax, int rhoCDmax, int rhoAB_CDmax,
       int rhoACmax, int rhoBDmax, int rhoAC_BDmax,
       u34utils::SO3_CGC& cgcs, u34utils::SO3_Racah& racahs,
       u3::WigU3U2_db& wigU3U2db)
  {
    std::vector<double> U(rhoACmax * rhoBDmax * rhoAB_CDmax * rhoABmax *
                              rhoCDmax * rhoAC_BDmax, 0.0);

    // couple u3AB and u3C
    auto [u3AB_C, rhomaxsAB_C] = u3::decomposition(u3AB, u3C);

    // couple u3AC and u3B
    auto [u3AC_B, rhomaxsAC_B] = u3::decomposition(u3AC, u3B);

    std::vector<u3::irrep> u3ABC;
    std::vector<int> rhoAB_Cmax, rhoAC_Bmax;

    for (int i = 0; i < u3AB_C.size(); i++)
    {
      for (int j = 0; j < u3AC_B.size(); j++)
      {
        if (u3AB_C[i] == u3AC_B[j])
        {
          u3ABC.push_back(u3AB_C[i]);
          rhoAB_Cmax.push_back(rhomaxsAB_C[i]);
          rhoAC_Bmax.push_back(rhomaxsAC_B[j]);
          break;
        }
      }
    }

    for (int i = 0; i < u3ABC.size(); i++)
    {
      std::vector<u3::gelfand> gArrABC = u3::gelfandBasis(u3ABC[i]);
/*
      auto [u3ABC_D, rhomaxsABC_D] = u3::decomposition(u3ABC[i], u3D);
      int rhoABC_Dmax = 0;
      for (int j = 0; j < u3ABC_D.size(); j++)
      {
        if (u3ABC_D[j] == u3ABCD)
        {
          rhoABC_Dmax = rhomaxsABC_D[j];
          break;
        }
      }
*/
      int rhoABC_Dmax = u3::Rho(u3ABC[i], u3D, u3ABCD);
      if (rhoABC_Dmax == 0)
        continue;

      //std::cout << u3ABC[i] << rhoAB_Cmax[i] << "," << rhoAC_Bmax[i] << ","
      // << rhoABC_Dmax << std::endl;

      std::vector<double> U6_1 =
          u3::U6(u3AC, gArrAC, u3B, gArrB, u3ABCD, gArrABCD, u3D, gArrD,
                 u3ABC[i], gArrABC, rhoAC_Bmax[i], rhoABC_Dmax, u3BD, gArrBD,
                 rhoBDmax, rhoAC_BDmax, cgcs, racahs, wigU3U2db);

      std::vector<double> U6_2 =
          u3::U6(u3AB, gArrAB, u3C, gArrC, u3ABCD, gArrABCD, u3D, gArrD,
                 u3ABC[i], gArrABC, rhoAB_Cmax[i], rhoABC_Dmax, u3CD, gArrCD,
                 rhoCDmax, rhoAB_CDmax, cgcs, racahs, wigU3U2db);

      std::vector<double> Z6 =
          u3::Z6(u3B, gArrB, u3A, gArrA, u3ABC[i], gArrABC, u3C, gArrC, u3AB,
                 gArrAB, rhoABmax, rhoAB_Cmax[i], u3AC, gArrAC, rhoACmax,
                 rhoAC_Bmax[i], cgcs, racahs, wigU3U2db);

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
                for (int rAC_BD = 0; rAC_BD < rhoAC_BDmax; rAC_BD++)
                {
                  for (int rAC_B = 0; rAC_B < rhoAC_Bmax[i]; rAC_B++)
                  {
                    for (int rAB_C = 0; rAB_C < rhoAB_Cmax[i]; rAB_C++)
                    {
                      for (int rABC_D = 0; rABC_D < rhoABC_Dmax; rABC_D++)
                      {
                        U[((((rAC * rhoBDmax + rBD) * rhoAB_CDmax + rAB_CD) *
                        rhoABmax + rAB) * rhoCDmax + rCD) * rhoAC_BDmax + rAC_BD] 
                        += U6_1[((rAC_B * rhoABC_Dmax + rABC_D) * rhoBDmax +
                                  rBD) * rhoAC_BDmax + rAC_BD] 
                         * Z6[((rAB * rhoAB_Cmax[i] + rAB_C) * rhoACmax +
                                rAC) * rhoAC_Bmax[i] + rAC_B] 
                         * U6_2[((rAB_C * rhoABC_Dmax + rABC_D) * rhoCDmax +
                                  rCD) * rhoAB_CDmax + rAB_CD];
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

    for (int i = 0; i < U.size(); i++)
      if (std::fabs(U[i]) < 1e-14)
        U[i] = 0.0;

    return U;
  }

} // namespace u3

#endif