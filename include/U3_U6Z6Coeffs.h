#ifndef U3LIB_U3_U6Z6COEFFS_H
#define U3LIB_U3_U6Z6COEFFS_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "U3basis.h"
#include "U3U2xU3U2.h"

#include <so3Coeffs/SO3_CGC.h>
#include <so3Coeffs/SO3_Racah.h>

namespace u3
{
  // input by users: u3A, u3B, u3ABC, u3C, u3AB, u3BC, all normalized
  inline std::vector<double> 
    U6(const u3::irrep& u3A, const std::vector<u3::gelfand>& gArrA,
       const u3::irrep& u3B, const std::vector<u3::gelfand>& gArrB,
       const u3::irrep& u3ABC, const std::vector<u3::gelfand>& gArrABC,
       const u3::irrep& u3C, const std::vector<u3::gelfand>& gArrC,
       const u3::irrep& u3AB, const std::vector<u3::gelfand>& gArrAB,
       int rhoABmax, int rhoAB_Cmax, 
       const u3::irrep& u3BC, const std::vector<u3::gelfand>& gArrBC, 
       int rhoBCmax, int rhoA_BCmax, 
       u34utils::SO3_CGC& cgcs, u34utils::SO3_Racah& racahs)
  {
    // coupling <gA;hwBC|hwABC> rhoA_BC
    auto coeffsA_BC = u3::WigU3U2(u3A, u3BC, u3ABC, rhoA_BCmax, 
      u3ABC.n13() + u3ABC.n23(), gArrA, gArrBC, gArrABC, cgcs); 

    // coupling <gAB;gC|hwABC> rhoAB_C
    auto coeffsAB_C = u3::WigU3U2(u3AB, u3C, u3ABC, rhoAB_Cmax, 
      u3ABC.n13() + u3ABC.n23(), gArrAB, gArrC, gArrABC, cgcs); 

    // coupling <gB;gC|hwBC> rhoBC
    auto coeffsBC = u3::WigU3U2(u3B, u3C, u3BC, rhoBCmax, 
      u3BC.n13() + u3BC.n23(), gArrB, gArrC, gArrBC, cgcs); 

    // coupling <gA;gB|gAB> rhoAB
    auto coeffsAB = u3::WigU3U2(u3A, u3B, u3AB, rhoABmax, -1, 
      gArrA, gArrB, gArrAB, cgcs); 

    Eigen::MatrixXd A(rhoA_BCmax, rhoA_BCmax);
    Eigen::MatrixXd B(rhoA_BCmax, rhoABmax * rhoAB_Cmax * rhoBCmax);
    B.setZero();

    int r = 0;

    for (int n12A = u3A.n13(); n12A >= u3A.n23(); n12A--)
    {
      for (int n22A = u3A.n23(); n22A >= u3A.n33(); n22A--)
      {
        auto accA_BC = coeffsA_BC.access(n12A, n22A, u3BC.n13(), u3BC.n23(),
                                         u3ABC.n13(), u3ABC.n23());
        if (accA_BC.allowed == true)
        {
          A.row(r) = coeffsA_BC.get_row(accA_BC);

          for (int n12B = u3B.n13(); n12B >= u3B.n23(); n12B--)
          {
            for (int n22B = u3B.n23(); n22B >= u3B.n33(); n22B--)
            {
              for (int n12AB = u3AB.n13(); n12AB >= u3AB.n23(); n12AB--)
              {
                for (int n22AB = u3AB.n23(); n22AB >= u3AB.n33(); n22AB--)
                {
                  auto accAB =
                      coeffsAB.access(n12A, n22A, n12B, n22B, n12AB, n22AB);
                  if (accAB.allowed == true)
                  {
                    for (int n12C = u3C.n13(); n12C >= u3C.n23(); n12C--)
                    {
                      for (int n22C = u3C.n23(); n22C >= u3C.n33(); n22C--)
                      {
                        auto accAB_C = coeffsAB_C.access(
                            n12AB, n22AB, n12C, n22C, u3ABC.n13(), u3ABC.n23());
                        auto accBC = coeffsBC.access(n12B, n22B, n12C, n22C,
                                                     u3BC.n13(), u3BC.n23());
                        if (accAB_C.allowed == true && accBC.allowed == true)
                        {
                          int llA = n12A - n22A, llB = n12B - n22B,
                              ll = u3ABC.n13() - u3ABC.n23(), llC = n12C - n22C,
                              llAB = n12AB - n22AB,
                              llBC = u3BC.n13() - u3BC.n23();

                          double racah =
                              racahs.compute(llA, llB, ll, llC, llAB, llBC);

                          for (int rAB = 0; rAB < rhoABmax; rAB++)
                            for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
                              for (int rBC = 0; rBC < rhoBCmax; rBC++)
                                B(r, (rAB * rhoAB_Cmax + rAB_C) * rhoBCmax +
                                         rBC) +=
                                    racah * coeffsAB.get(accAB, rAB) *
                                    coeffsAB_C.get(accAB_C, rAB_C) *
                                    coeffsBC.get(accBC, rBC);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          r++;
        }
      }
      if (r == rhoA_BCmax)
        break;
    }

    // std::cout << A << std::endl;
    // std::cout << B << std::endl;

    Eigen::MatrixXd X = A.householderQr().solve(B);

    std::vector<double> U(rhoABmax * rhoAB_Cmax * rhoBCmax * rhoA_BCmax);

    for (int rAB = 0; rAB < rhoABmax; rAB++)
      for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
        for (int rBC = 0; rBC < rhoBCmax; rBC++)
          for (int rA_BC = 0; rA_BC < rhoA_BCmax; rA_BC++)
          {
            if (std::fabs(X(rA_BC, (rAB * rhoAB_Cmax + rAB_C) * rhoBCmax +
                                       rBC)) < 1e-14)
              X(rA_BC, (rAB * rhoAB_Cmax + rAB_C) * rhoBCmax + rBC) = 0.0;

            U[((rAB * rhoAB_Cmax + rAB_C) * rhoBCmax + rBC) * rhoA_BCmax +
              rA_BC] = X(rA_BC, (rAB * rhoAB_Cmax + rAB_C) * rhoBCmax + rBC);
          }

    return U;
  }

  // input by users: u3B, u3B, u3ABC, u3C, u3AB, u3AC
  inline std::vector<double> 
    Z6(const u3::irrep& u3B, const std::vector<u3::gelfand>& gArrB,
       const u3::irrep& u3A, const std::vector<u3::gelfand>& gArrA,
       const u3::irrep& u3ABC, const std::vector<u3::gelfand>& gArrABC,
       const u3::irrep& u3C, const std::vector<u3::gelfand>& gArrC,
       const u3::irrep& u3AB, const std::vector<u3::gelfand>& gArrAB,
       int rhoABmax, int rhoAB_Cmax, 
       const u3::irrep& u3AC, const std::vector<u3::gelfand>& gArrAC, 
       int rhoACmax, int rhoAC_Bmax, 
       u34utils::SO3_CGC& cgcs, u34utils::SO3_Racah& racahs)
  {
    // coupling <hwAC;gB|hwABC> rhoAC_B
    auto coeffsAC_B = u3::WigU3U2(u3AC, u3B, u3ABC, rhoAC_Bmax, 
      u3ABC.n13() + u3ABC.n23(), gArrAC, gArrB, gArrABC, cgcs);

    // coupling <gAB;gC|hwABC> rhoAB_C
    auto coeffsAB_C = u3::WigU3U2(u3AB, u3C, u3ABC, rhoAB_Cmax, 
      u3ABC.n13() + u3ABC.n23(), gArrAB, gArrC, gArrABC, cgcs);

    // coupling <gA;gC|hwAC> rhoAC
    auto coeffsAC = u3::WigU3U2(u3A, u3C, u3AC, rhoACmax, 
      u3AC.n13() + u3AC.n23(), gArrA, gArrC, gArrAC, cgcs);

    // coupling <gA;gB|gAB> rhoAB
    auto coeffsAB = u3::WigU3U2(u3A, u3B, u3AB, rhoABmax, -1, 
      gArrA, gArrB, gArrAB, cgcs);

    Eigen::MatrixXd A(rhoAC_Bmax, rhoAC_Bmax);
    Eigen::MatrixXd B(rhoAC_Bmax, rhoABmax * rhoAB_Cmax * rhoACmax);
    B.setZero();

    int r = 0;

    //std::cout << "u3AC: " << u3AC << ", u3B: " << u3B << ", u3ABC: " << u3ABC << std::endl;
    for (int n12B = u3B.n13(); n12B >= u3B.n23(); n12B--)
    {
      for (int n22B = u3B.n23(); n22B >= u3B.n33(); n22B--)
      {
        auto accAC_B = coeffsAC_B.access(u3AC.n13(), u3AC.n23(), n12B, n22B,
                                         u3ABC.n13(), u3ABC.n23());
        if (accAC_B.allowed == true)
        {
          A.row(r) = coeffsAC_B.get_row(accAC_B);

          for (int n12A = u3A.n13(); n12A >= u3A.n23(); n12A--)
          {
            for (int n22A = u3A.n23(); n22A >= u3A.n33(); n22A--)
            {
              for (int n12AB = u3AB.n13(); n12AB >= u3AB.n23(); n12AB--)
              {
                for (int n22AB = u3AB.n23(); n22AB >= u3AB.n33(); n22AB--)
                {
                  auto accAB =
                      coeffsAB.access(n12A, n22A, n12B, n22B, n12AB, n22AB);
                  if (accAB.allowed == true)
                  {
                    for (int n12C = u3C.n13(); n12C >= u3C.n23(); n12C--)
                    {
                      for (int n22C = u3C.n23(); n22C >= u3C.n33(); n22C--)
                      {
                        auto accAB_C = coeffsAB_C.access(
                            n12AB, n22AB, n12C, n22C, u3ABC.n13(), u3ABC.n23());
                        auto accAC = coeffsAC.access(n12A, n22A, n12C, n22C,
                                                     u3AC.n13(), u3AC.n23());
                        if (accAB_C.allowed == true && accAC.allowed == true)
                        {
                          int llA = n12A - n22A, llB = n12B - n22B,
                              ll = u3ABC.n13() - u3ABC.n23(), llC = n12C - n22C,
                              llAB = n12AB - n22AB,
                              llAC = u3AC.n13() - u3AC.n23();

                          double phase =
                              ((llA + ll - llAB - llAC) % 4 == 0) ? 1 : -1;
                          double racah =
                              racahs.compute(llB, llA, ll, llC, llAB, llAC);

                          for (int rAB = 0; rAB < rhoABmax; rAB++)
                            for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
                              for (int rAC = 0; rAC < rhoACmax; rAC++)
                                B(r, (rAB * rhoAB_Cmax + rAB_C) * rhoACmax +
                                         rAC) +=
                                    phase * racah * coeffsAB.get(accAB, rAB) *
                                    coeffsAB_C.get(accAB_C, rAB_C) *
                                    coeffsAC.get(accAC, rAC);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          r++;
        }
      }
      if (r == rhoAC_Bmax)
        break;
    }

    // std::cout << A << std::endl;
    // std::cout << B << std::endl;

    Eigen::MatrixXd X = A.householderQr().solve(B);

    std::vector<double> Z(rhoABmax * rhoAB_Cmax * rhoACmax * rhoAC_Bmax);

    for (int rAB = 0; rAB < rhoABmax; rAB++)
      for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
        for (int rAC = 0; rAC < rhoACmax; rAC++)
          for (int rAC_B = 0; rAC_B < rhoAC_Bmax; rAC_B++)
          {
            if (std::fabs(X(rAC_B, (rAB * rhoAB_Cmax + rAB_C) * rhoACmax +
                                       rAC)) < 1e-14)
              X(rAC_B, (rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC) = 0.0;

            Z[((rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC) * rhoAC_Bmax +
              rAC_B] = X(rAC_B, (rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC);
          }

    return Z;
  }

} // namespace u3
#endif