#ifndef U3LIB_U3SO3xU3SO3_SPECIFIC_H
#define U3LIB_U3SO3xU3SO3_SPECIFIC_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "cgc_so3/CGC_SO3.h"
#include "U3basis.h"
#include "U3SO3_to_U3U2U1.h"
#include "U3U2xU3U2.h"

namespace u3
{
  // generate Wigner U3->SO(3) coefficients for specific l1, l2, l3
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  WigU3SO3_specific(int& rhomax, u3::physicalvector& pA,
          u3::physicalvector& pB, u3::physicalvector& pC)
  {

    // find all the reduced canonical coupling coefficients
    std::vector<u3::gelfand> gArrA, gArrB, gArrC;

    auto canoCoeffs =
      u3::WigU3U2(pA.toprow, pB.toprow, pC.toprow, rhomax, gArrA, gArrB, gArrC);
    // note that u3C and all gelfand states are already normalized here

    // find transformation for the first state
    int mA = pC.L - pB.L;
    std::vector<int> canoBasisA;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
      transformationA = u3::phys_to_cano(pA, mA, canoBasisA, gArrA);

    // find transformation for the second state
    std::vector<int> canoBasisB;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
      transformationB = u3::phys_to_cano(pB, pB.L, canoBasisB, gArrB);

    // find transformation for the third state
    std::vector<int> canoBasisC;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> 
      transformationC = u3::phys_to_cano(pC, pC.L, canoBasisC, gArrC);

    u3::CGC_SO3 cgcs;
    double cgc1 = cgcs.compute(2*pA.L, 2*pB.L, 2*pC.L, 2*mA, 2*pB.L, 2*pC.L);

    Eigen::MatrixXd coeffs(pA.kappa * pB.kappa * pC.kappa, rhomax);
    coeffs.setZero();

    for (int iA = 0; iA < canoBasisA.size(); iA++)
    {
      for (int iB = 0; iB < canoBasisB.size(); iB++)
      {
        for (int iC = 0; iC < canoBasisC.size(); iC++)
        {
          Eigen::MatrixXd wigU3U2U1;
          double cgc2;

          int n12A = gArrA[canoBasisA[iA]].n12(), n22A = gArrA[canoBasisA[iA]].n22(),
              n12B = gArrB[canoBasisB[iB]].n12(), n22B = gArrB[canoBasisB[iB]].n22(),
              n12C = gArrC[canoBasisC[iC]].n12(), n22C = gArrC[canoBasisC[iC]].n22();

          auto acc = canoCoeffs.access(n12A, n22A, n12B, n22B, n12C, n22C);

          int jj1 = gArrA[canoBasisA[iA]].Lambdax2(), mm1 = gArrA[canoBasisA[iA]].MLambdax2(),
              jj2 = gArrB[canoBasisB[iB]].Lambdax2(), mm2 = gArrB[canoBasisB[iB]].MLambdax2(),
              jj3 = gArrC[canoBasisC[iC]].Lambdax2(), mm3 = gArrC[canoBasisC[iC]].MLambdax2();

          if ((acc.allowed == false) || ((mm1 + mm2) != mm3))
            continue;

          else
          {
            wigU3U2U1 = canoCoeffs.get_row(acc);
            cgc2 = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);
          }

          for (int r = 0; r < rhomax; r++)
            for (int k1 = 0; k1 < pA.kappa; k1++)
              for (int k2 = 0; k2 < pB.kappa; k2++)
                for (int k3 = 0; k3 < pC.kappa; k3++)
                  coeffs((k1 * pB.kappa + k2) * pC.kappa + k3, r) +=
                    transformationA(iA,k1) * transformationB(iB,k2) *
                    transformationC(iC,k3) * wigU3U2U1(0,r) * cgc2;
        }
      }
    }

    return coeffs/cgc1;
  }

} // namespace u3

#endif
