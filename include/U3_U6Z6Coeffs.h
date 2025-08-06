#ifndef U3LIB_U3_U6Z6COEFFS_H
#define U3LIB_U3_U6Z6COEFFS_H

#include <cmath>
#include <iostream>
#include <map>
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
     int rhoABmax, int rhoAB_Cmax, const u3::irrep& u3BC,
     const std::vector<u3::gelfand>& gArrBC, int rhoBCmax, int rhoA_BCmax,
     qcx_utils::SO3_CGC& cgcs, qcx_utils::SO3_Racah& racahs,
     u3::WigU3U2_db& wigU3U2db)
  {
    Eigen::MatrixXd A(rhoA_BCmax, rhoA_BCmax);
    Eigen::MatrixXd B(rhoA_BCmax, rhoABmax * rhoAB_Cmax * rhoBCmax);
    B.setZero();

    int r = 0;
    bool done = false;
    for (int n12A = u3A.n13(); n12A >= u3A.n23() && !done; n12A--)
    {
      for (int n22A = u3A.n23(); n22A >= u3A.n33() && !done; n22A--)
      {
        auto [allowA_BC, accA_BC] = wigU3U2db.get(
            u3A.n13(), u3A.n23(), u3A.n33(), n12A, n22A,
            u3BC.n13(), u3BC.n23(), u3BC.n33(), u3BC.n13(), u3BC.n23(),
            u3ABC.n13(), u3ABC.n23(), u3ABC.n33(), u3ABC.n13(), u3ABC.n23());

        if (allowA_BC == true)
        {
          A.row(r) = accA_BC;

          for (int n12B = u3B.n13(); n12B >= u3B.n23(); n12B--)
          {
            for (int n22B = u3B.n23(); n22B >= u3B.n33(); n22B--)
            {
              for (int n12AB = u3AB.n13(); n12AB >= u3AB.n23(); n12AB--)
              {
                for (int n22AB = u3AB.n23(); n22AB >= u3AB.n33(); n22AB--)
                {
                  auto [allowAB, accAB] = wigU3U2db.get(
                      u3A.n13(), u3A.n23(), u3A.n33(), n12A, n22A,
                      u3B.n13(), u3B.n23(), u3B.n33(), n12B, n22B,
                      u3AB.n13(), u3AB.n23(), u3AB.n33(), n12AB, n22AB);
                  
                  if (allowAB == true)
                  {
                    for (int n12C = u3C.n13(); n12C >= u3C.n23(); n12C--)
                    {
                      for (int n22C = u3C.n23(); n22C >= u3C.n33(); n22C--)
                      {
                        auto [allowAB_C, accAB_C] = wigU3U2db.get(
                          u3AB.n13(), u3AB.n23(), u3AB.n33(), n12AB, n22AB,
                          u3C.n13(), u3C.n23(), u3C.n33(), n12C, n22C,
                          u3ABC.n13(), u3ABC.n23(), u3ABC.n33(), u3ABC.n13(), u3ABC.n23());

                        auto [allowBC, accBC] = wigU3U2db.get(
                          u3B.n13(), u3B.n23(), u3B.n33(), n12B, n22B,
                          u3C.n13(), u3C.n23(), u3C.n33(), n12C, n22C,
                          u3BC.n13(), u3BC.n23(), u3BC.n33(), u3BC.n13(), u3BC.n23());

                        if (allowAB_C == true && allowBC == true)
                        {
                          int llA = n12A - n22A, ll = u3ABC.n13() - u3ABC.n23(), 
                              llBC = u3BC.n13() - u3BC.n23(), llB = n12B - n22B, 
                              llAB = n12AB - n22AB, llC = n12C - n22C;
                          double racah =
                              racahs.compute(llA, llB, ll, llC, llAB, llBC);

                          for (int rAB = 0; rAB < rhoABmax; rAB++)
                            for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
                              for (int rBC = 0; rBC < rhoBCmax; rBC++)
                                B(r, (rAB * rhoAB_Cmax + rAB_C) * rhoBCmax + rBC) 
                                  += racah * accAB(rAB) * accAB_C(rAB_C) * accBC(rBC);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          r++;

          if (r == rhoA_BCmax)
            done = true;
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
            int col = (rAB * rhoAB_Cmax + rAB_C) * rhoBCmax + rBC;
            if (std::fabs(X(rA_BC, col)) < 1e-14)
              X(rA_BC, col) = 0.0;

            U[((rAB * rhoAB_Cmax + rAB_C) * rhoBCmax + rBC) * rhoA_BCmax +
              rA_BC] = X(rA_BC, col);
          }

    return U;
  }

#define U3LIB_ENABLE_U3U6_CACHE
  class U3U6_db
  {
#ifdef U3LIB_ENABLE_U3U6_CACHE
    using key_type = std::array<int, 12>;
    std::map<key_type, std::vector<double>> map_;
#endif

    qcx_utils::SO3_CGC& cgcs_;
    qcx_utils::SO3_Racah& racahs_;
    u3::WigU3U2_db& wigU3U2db_;

    auto calculate(const u3::irrep& u3A, const u3::irrep& u3B, 
                   const u3::irrep& u3ABC, const u3::irrep& u3C,
                   const u3::irrep& u3AB, int rhoABmax, int rhoAB_Cmax,
                   const u3::irrep& u3BC, int rhoBCmax, int rhoA_BCmax)
    {
      std::vector<u3::gelfand> gArrA = u3::gelfandBasis(u3A),
                               gArrB = u3::gelfandBasis(u3B),
                               gArrC = u3::gelfandBasis(u3C),
                               gArrAB = u3::gelfandBasis(u3AB),
                               gArrBC = u3::gelfandBasis(u3BC),
                               gArrABC = u3::gelfandBasis(u3ABC);

      return U6(u3A,gArrA,u3B,gArrB,u3ABC,gArrABC,u3C,gArrC,
                u3AB,gArrAB,rhoABmax,rhoAB_Cmax,
                u3BC,gArrBC,rhoBCmax,rhoA_BCmax,
                cgcs_,racahs_,wigU3U2db_);
    }

  public:
    U3U6_db(qcx_utils::SO3_CGC& cgcs, qcx_utils::SO3_Racah& racahs,
                                          u3::WigU3U2_db& wigU3U2db) 
      : cgcs_(cgcs), racahs_(racahs), wigU3U2db_(wigU3U2db)
    {}

    auto get(int n13a, int n23a, int n33a, int n13b, int n23b, int n33b, 
      int n13abc, int n23abc, int n33abc, int n13c, int n23c, int n33c, 
      int n13ab, int n23ab, int n33ab, int n13bc, int n23bc, int n33bc)
    {
      u3::irrep u3A(n13a - n33a, n23a - n33a, 0), 
                u3B(n13b - n33b, n23b - n33b, 0), 
                u3C(n13c - n33c, n23c - n33c, 0), 
                u3AB(n13ab - n33ab, n23ab - n33ab, 0), 
                u3BC(n13bc - n33bc, n23bc - n33bc, 0), 
                u3ABC(n13abc - n33abc, n23abc - n33abc, 0);

      int rhoABmax = Rho(u3A,u3B,u3AB), rhoAB_Cmax = Rho(u3AB,u3C,u3ABC),
          rhoBCmax = Rho(u3B,u3C,u3BC), rhoA_BCmax = Rho(u3A,u3BC,u3ABC);

      bool flag = false;
      std::vector<double> c;

      if (rhoABmax > 0 && rhoAB_Cmax > 0 && rhoBCmax > 0 && rhoA_BCmax > 0)
      {
#ifdef U3LIB_ENABLE_U3U6_CACHE
        key_type key = {u3A.n13(), u3A.n23(),
                        u3B.n13(), u3B.n23(),
                        u3ABC.n13(), u3ABC.n23(),
                        u3C.n13(), u3C.n23(),
                        u3AB.n13(), u3AB.n23(),
                        u3BC.n13(), u3BC.n23()};
        auto iter = map_.find(key);
        if (iter == map_.end())
        {
          auto p = map_.try_emplace(key, calculate(u3A, u3B, u3ABC, u3C, 
                    u3AB, rhoABmax, rhoAB_Cmax, u3BC, rhoBCmax, rhoA_BCmax));
          assert(p.second);
          iter = p.first;
        }

        c = iter->second;
#else
        c = calculate(u3A, u3B, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                            u3BC, rhoBCmax, rhoA_BCmax);
#endif
        flag = true;
      }

      return std::make_pair(flag, c);
    }
  };

  // input by users: u3B, u3A, u3ABC, u3C, u3AB, u3AC
  inline std::vector<double>
  Z6(const u3::irrep& u3B, const std::vector<u3::gelfand>& gArrB,
     const u3::irrep& u3A, const std::vector<u3::gelfand>& gArrA,
     const u3::irrep& u3ABC, const std::vector<u3::gelfand>& gArrABC,
     const u3::irrep& u3C, const std::vector<u3::gelfand>& gArrC,
     const u3::irrep& u3AB, const std::vector<u3::gelfand>& gArrAB,
     int rhoABmax, int rhoAB_Cmax, const u3::irrep& u3AC,
     const std::vector<u3::gelfand>& gArrAC, int rhoACmax, int rhoAC_Bmax,
     qcx_utils::SO3_CGC& cgcs, qcx_utils::SO3_Racah& racahs,
     u3::WigU3U2_db& wigU3U2db)
  {
    Eigen::MatrixXd A(rhoAC_Bmax, rhoAC_Bmax);
    Eigen::MatrixXd B(rhoAC_Bmax, rhoABmax * rhoAB_Cmax * rhoACmax);
    B.setZero();

    int r = 0;
    bool done = false;
    for (int n12B = u3B.n13(); n12B >= u3B.n23() && !done; n12B--)
    {
      for (int n22B = u3B.n23(); n22B >= u3B.n33() && !done; n22B--)
      {
        auto [allowAC_B, accAC_B] = wigU3U2db.get(
          u3AC.n13(), u3AC.n23(), u3AC.n33(), u3AC.n13(), u3AC.n23(), 
          u3B.n13(), u3B.n23(), u3B.n33(), n12B, n22B,
          u3ABC.n13(), u3ABC.n23(), u3ABC.n33(), u3ABC.n13(), u3ABC.n23());

        if (allowAC_B == true)
        {
          A.row(r) = accAC_B;

          for (int n12A = u3A.n13(); n12A >= u3A.n23(); n12A--)
          {
            for (int n22A = u3A.n23(); n22A >= u3A.n33(); n22A--)
            {
              for (int n12AB = u3AB.n13(); n12AB >= u3AB.n23(); n12AB--)
              {
                for (int n22AB = u3AB.n23(); n22AB >= u3AB.n33(); n22AB--)
                {
                  auto [allowAB, accAB] = wigU3U2db.get(
                    u3A.n13(), u3A.n23(), u3A.n33(), n12A, n22A, 
                    u3B.n13(), u3B.n23(), u3B.n33(), n12B, n22B, 
                    u3AB.n13(), u3AB.n23(), u3AB.n33(), n12AB, n22AB);

                  if (allowAB == true)
                  {
                    for (int n12C = u3C.n13(); n12C >= u3C.n23(); n12C--)
                    {
                      for (int n22C = u3C.n23(); n22C >= u3C.n33(); n22C--)
                      {
                        auto [allowAB_C, accAB_C] = wigU3U2db.get(
                            u3AB.n13(), u3AB.n23(), u3AB.n33(), n12AB, n22AB, 
                            u3C.n13(), u3C.n23(), u3C.n33(), n12C, n22C, 
                            u3ABC.n13(), u3ABC.n23(), u3ABC.n33(), u3ABC.n13(), u3ABC.n23());
                        
                        auto [allowAC, accAC] = wigU3U2db.get(
                          u3A.n13(), u3A.n23(), u3A.n33(), n12A, n22A, 
                          u3C.n13(), u3C.n23(), u3C.n33(), n12C, n22C,
                          u3AC.n13(), u3AC.n23(), u3AC.n33(), u3AC.n13(), u3AC.n23());

                        if (allowAB_C == true && allowAC == true)
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
                                B(r, (rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC) 
                                  += phase * racah * accAB(rAB) * accAB_C(rAB_C) 
                                  * accAC(rAC);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          r++;

          if (r == rhoAC_Bmax)
            done = true;
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
            int col = (rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC;
            if (std::fabs(X(rAC_B, col)) < 1e-14)
              X(rAC_B, col) = 0.0;

            Z[((rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC) * rhoAC_Bmax +
              rAC_B] = X(rAC_B, col);
          }

    return Z;
  }

#define U3LIB_ENABLE_U3Z6_CACHE
  class U3Z6_db
  {
#ifdef U3LIB_ENABLE_U3Z6_CACHE
    using key_type = std::array<int, 12>;
    std::map<key_type, std::vector<double>> map_;
#endif

    qcx_utils::SO3_CGC& cgcs_;
    qcx_utils::SO3_Racah& racahs_;
    u3::WigU3U2_db& wigU3U2db_;

    auto calculate(const u3::irrep& u3B, const u3::irrep& u3A, 
                   const u3::irrep& u3ABC, const u3::irrep& u3C,
                   const u3::irrep& u3AB, int rhoABmax, int rhoAB_Cmax,
                   const u3::irrep& u3AC, int rhoACmax, int rhoAC_Bmax)
    {
      std::vector<u3::gelfand> gArrA = u3::gelfandBasis(u3A),
                               gArrB = u3::gelfandBasis(u3B),
                               gArrC = u3::gelfandBasis(u3C),
                               gArrAB = u3::gelfandBasis(u3AB),
                               gArrAC = u3::gelfandBasis(u3AC),
                               gArrABC = u3::gelfandBasis(u3ABC);

      return U6(u3B,gArrB,u3A,gArrA,u3ABC,gArrABC,u3C,gArrC,
                u3AB,gArrAB,rhoABmax,rhoAB_Cmax,
                u3AC,gArrAC,rhoACmax,rhoAC_Bmax,
                cgcs_,racahs_,wigU3U2db_);
    }

  public:
    U3Z6_db(qcx_utils::SO3_CGC& cgcs, qcx_utils::SO3_Racah& racahs,
                                          u3::WigU3U2_db& wigU3U2db) 
      : cgcs_(cgcs), racahs_(racahs), wigU3U2db_(wigU3U2db)
    {}

    auto get(int n13b, int n23b, int n33b, int n13a, int n23a, int n33a,
      int n13abc, int n23abc, int n33abc, int n13c, int n23c, int n33c, 
      int n13ab, int n23ab, int n33ab, int n13ac, int n23ac, int n33ac)
    {
      u3::irrep u3A(n13a - n33a, n23a - n33a, 0), 
                u3B(n13b - n33b, n23b - n33b, 0), 
                u3C(n13c - n33c, n23c - n33c, 0), 
                u3AB(n13ab - n33ab, n23ab - n33ab, 0), 
                u3AC(n13ac - n33ac, n23ac - n33ac, 0), 
                u3ABC(n13abc - n33abc, n23abc - n33abc, 0);

      int rhoABmax = Rho(u3A,u3B,u3AB), rhoAB_Cmax = Rho(u3AB,u3C,u3ABC),
          rhoACmax = Rho(u3A,u3C,u3AC), rhoAC_Bmax = Rho(u3AC,u3B,u3ABC);

      bool flag = false;
      std::vector<double> c;

      if (rhoABmax > 0 && rhoAB_Cmax > 0 && rhoACmax > 0 && rhoAC_Bmax > 0)
      {
#ifdef U3LIB_ENABLE_U3U6_CACHE
        key_type key = {u3B.n13(), u3B.n23(),
                        u3A.n13(), u3A.n23(),
                        u3ABC.n13(), u3ABC.n23(),
                        u3C.n13(), u3C.n23(),
                        u3AB.n13(), u3AB.n23(),
                        u3AC.n13(), u3AC.n23()};
        auto iter = map_.find(key);
        if (iter == map_.end())
        {
          auto p = map_.try_emplace(key, calculate(u3B, u3A, u3ABC, u3C, 
                    u3AB, rhoABmax, rhoAB_Cmax, u3AC, rhoACmax, rhoAC_Bmax));
          assert(p.second);
          iter = p.first;
        }

        c = iter->second;
#else
        c = calculate(u3B, u3A, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                            u3AC, rhoACmax, rhoAC_Bmax);
#endif
        flag = true;
      }

      return std::make_pair(flag, c);
    }
  };

} // namespace u3
#endif
