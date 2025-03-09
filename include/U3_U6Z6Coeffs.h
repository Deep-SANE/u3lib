#ifndef U3LIB_U3_U6Z6COEFFS_H
#define U3LIB_U3_U6Z6COEFFS_H

#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>

#include "so3Coeffs/SO3_CGC.h"
#include "so3Coeffs/SO3_Racah.h"
#include "U3basis.h"
#include "U3U2xU3U2.h"

namespace u3
{
  // input by users: u3A, u3B, u3ABC, u3C, u3AB, u3BC
  std::vector<double> U6(u3::irrep& u3A, std::vector<u3::gelfand>& gArrA,
                         u3::irrep& u3B, std::vector<u3::gelfand>& gArrB,
                         u3::irrep& u3ABC, std::vector<u3::gelfand>& gArrABC,
                         u3::irrep& u3C, std::vector<u3::gelfand>& gArrC,
                         u3::irrep& u3AB, std::vector<u3::gelfand>& gArrAB,
                         int rhoABmax, int rhoAB_Cmax, u3::irrep& u3BC,
                         std::vector<u3::gelfand>& gArrBC, int rhoBCmax,
                         int rhoA_BCmax, u3::SO3_CGC& cgcs,
                         u3::SO3_Racah& racahs)
  {
    // coupling <gA;hwBC|hwABC> rhoA_BC
    u3::WigU3U2_coeffs coeffsA_BC(u3A, u3BC, u3ABC, rhoA_BCmax);
    u3::HighestWeightReducedCoeffs(gArrA, gArrBC, gArrABC, cgcs, coeffsA_BC);

    // coupling <gAB;gC|hwABC> rhoAB_C
    u3::WigU3U2_coeffs coeffsAB_C(u3AB, u3C, u3ABC, rhoAB_Cmax);
    u3::HighestWeightReducedCoeffs(gArrAB, gArrC, gArrABC, cgcs, coeffsAB_C);

    // coupling <gB;gC|hwBC> rhoBC
    u3::WigU3U2_coeffs coeffsBC(u3B, u3C, u3BC, rhoBCmax);
    u3::HighestWeightReducedCoeffs(gArrB, gArrC, gArrBC, cgcs, coeffsBC);

    // coupling <gA;gB|gAB> rhoAB
    u3::WigU3U2_coeffs coeffsAB(u3A, u3B, u3AB, rhoABmax);
    u3::HighestWeightReducedCoeffs(gArrA, gArrB, gArrAB, cgcs, coeffsAB);

    int w2max = u3AB.n13() + u3AB.n23(), w2min = u3AB.n23() + u3AB.n33();
    std::vector<std::vector<int>> all_same_weight_gArrAB(w2max - w2min + 1);
    for (int i = 0; i < gArrAB.size(); i++)
      if (gArrAB[i].n12() == gArrAB[i].n11())
        all_same_weight_gArrAB[w2max - gArrAB[i].n12() - gArrAB[i].n22()]
            .push_back(i);

    for (int w2 = w2max - 1; w2 >= w2min; w2--)
    {
      u3::LowerWeightReducedCoeffs(gArrA, gArrB, gArrAB, rhoABmax,
                                   all_same_weight_gArrAB[w2max - w2], cgcs,
                                   coeffsAB);
    }

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
  std::vector<double> Z6(u3::irrep& u3B, std::vector<u3::gelfand>& gArrB,
                         u3::irrep& u3A, std::vector<u3::gelfand>& gArrA,
                         u3::irrep& u3ABC, std::vector<u3::gelfand>& gArrABC,
                         u3::irrep& u3C, std::vector<u3::gelfand>& gArrC,
                         u3::irrep& u3AB, std::vector<u3::gelfand>& gArrAB,
                         int rhoABmax, int rhoAB_Cmax, u3::irrep& u3AC,
                         std::vector<u3::gelfand>& gArrAC, int rhoACmax,
                         int rhoAC_Bmax, u3::SO3_CGC& cgcs,
                         u3::SO3_Racah& racahs)
  {
    // coupling <hwAC;gB|hwABC> rhoAC_B
    u3::WigU3U2_coeffs coeffsAC_B(u3AC, u3B, u3ABC, rhoAC_Bmax);
    u3::HighestWeightReducedCoeffs(gArrAC, gArrB, gArrABC, cgcs, coeffsAC_B);

    // coupling <gAB;gC|hwABC> rhoAB_C
    u3::WigU3U2_coeffs coeffsAB_C(u3AB, u3C, u3ABC, rhoAB_Cmax);
    u3::HighestWeightReducedCoeffs(gArrAB, gArrC, gArrABC, cgcs, coeffsAB_C);

    // coupling <gA;gC|hwAC> rhoAC
    u3::WigU3U2_coeffs coeffsAC(u3A, u3C, u3AC, rhoACmax);
    u3::HighestWeightReducedCoeffs(gArrA, gArrC, gArrAC, cgcs, coeffsAC);

    // coupling <gA;gB|gAB> rhoAB
    u3::WigU3U2_coeffs coeffsAB(u3A, u3B, u3AB, rhoABmax);
    u3::HighestWeightReducedCoeffs(gArrA, gArrB, gArrAB, cgcs, coeffsAB);

    int w2max = u3AB.n13() + u3AB.n23(), w2min = u3AB.n23() + u3AB.n33();
    std::vector<std::vector<int>> all_same_weight_gArrAB(w2max - w2min + 1);
    for (int i = 0; i < gArrAB.size(); i++)
      if (gArrAB[i].n12() == gArrAB[i].n11())
        all_same_weight_gArrAB[w2max - gArrAB[i].n12() - gArrAB[i].n22()]
            .push_back(i);

    for (int w2 = w2max - 1; w2 >= w2min; w2--)
    {
      u3::LowerWeightReducedCoeffs(gArrA, gArrB, gArrAB, rhoABmax,
                                   all_same_weight_gArrAB[w2max - w2], cgcs,
                                   coeffsAB);
    }

    Eigen::MatrixXd A(rhoAC_Bmax, rhoAC_Bmax);
    Eigen::MatrixXd B(rhoAC_Bmax, rhoABmax * rhoAB_Cmax * rhoACmax);
    B.setZero();

    int r = 0;

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