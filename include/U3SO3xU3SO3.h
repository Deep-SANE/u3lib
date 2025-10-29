#ifndef U3LIB_U3SO3xU3SO3_H
#define U3LIB_U3SO3xU3SO3_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <tuple>
#include <vector>

#include <Eigen/Dense>

#include "U3basis.h"
#include "U3SO3_to_U3U2U1.h"
#include "U3U2xU3U2.h"

#include <so3Coeffs/SO3_CGC.h>

namespace u3
{
  // generate Wigner U3->SO(3) coefficients for multiple values of l1, l2, l3
  inline std::vector<Eigen::MatrixXd>
  WigU3SO3(int& rhomax, std::vector<u3::physicalvector>& pArrA,
           std::vector<u3::physicalvector>& pArrB,
           std::vector<u3::physicalvector>& pArrC,
           std::vector<std::tuple<int, int, int>>& tArr,
           qcx_utils::SO3_CGC& cgcs, u3::WigU3U2_db& wigU3U2db)
  {

    // build Gelfand basis
    std::vector<u3::gelfand> gArrA = u3::gelfandBasis(pArrA[0].toprow),
                             gArrB = u3::gelfandBasis(pArrB[0].toprow),
                             gArrC = u3::gelfandBasis(pArrC[0].toprow);

    // find transformation for physical states belong to the 2nd irrep
    std::vector<std::vector<int>> canoBasisB(pArrB.size());
    std::vector<Eigen::MatrixXd> transformationB(pArrB.size());
    for (int i = 0; i < pArrB.size(); i++)
      transformationB[i] =
          u3::phys_to_cano(pArrB[i], pArrB[i].L, canoBasisB[i], gArrB);

    // find transformation for physical states belong to the 3rd irrep
    std::vector<int> w2mins(pArrC.size());
    std::vector<std::vector<int>> canoBasisC(pArrC.size());
    std::vector<Eigen::MatrixXd> transformationC(pArrC.size());
    for (int i = 0; i < pArrC.size(); i++)
    {
      transformationC[i] =
          u3::phys_to_cano(pArrC[i], pArrC[i].L, canoBasisC[i], gArrC);
      w2mins[i] =
          gArrC[canoBasisC[i].back()].n12() + gArrC[canoBasisC[i].back()].n22();
    }

    // find all the reduced canonical coupling coefficients
    //int w2min = *std::min_element(w2mins.begin(), w2mins.end());
    //auto canoCoeffs =
    //    u3::WigU3U2(pArrA[0].toprow, pArrB[0].toprow, pArrC[0].toprow, rhomax,
    //                w2min, gArrA, gArrB, gArrC, cgcs);
    // note that u3C and all gelfand states are already normalized here

    std::vector<Eigen::MatrixXd> coeffs;

    int i, j, k;
    for (i = 0; i < pArrA.size(); i++)
    {
      std::vector<std::vector<int>> canoBasisA(1);
      std::vector<Eigen::MatrixXd> transformationA(1);
      transformationA[0] =
          u3::phys_to_cano(pArrA[i], pArrA[i].L, canoBasisA[0], gArrA);

      for (j = 0; j < pArrB.size(); j++)
      {
        for (k = 0; k < pArrC.size(); k++)
        {
          if (std::fabs(pArrA[i].L - pArrB[j].L) <= pArrC[k].L &&
              (pArrA[i].L + pArrB[j].L) >= pArrC[k].L)
          {
            tArr.emplace_back(i, j, k);

            Eigen::MatrixXd coeffs_temp(pArrA[i].kappa * pArrB[j].kappa 
                                                * pArrC[k].kappa, rhomax);
            coeffs_temp.setZero();

            int mA = pArrC[k].L - pArrB[j].L, mA_abs;
            if (mA >= 0)
              mA_abs = mA;
            else
              mA_abs = -mA;

            while ((pArrA[i].L - mA_abs + 1) > canoBasisA.size())
            {
              int m_new = pArrA[i].L - canoBasisA.size();
              std::vector<int> canoBasisA_new =
                  u3::expand_to_gelfand(m_new, gArrA);

              Eigen::MatrixXd transformationA_new = u3::lowerM(pArrA[i].kappa, 
                pArrA[i].L, m_new, gArrA, transformationA.back(), 
                canoBasisA.back(), canoBasisA_new);

              canoBasisA.push_back(canoBasisA_new);
              transformationA.push_back(transformationA_new);
            }

            double cgc1 =
                cgcs.compute(2 * pArrA[i].L, 2 * pArrB[j].L, 2 * pArrC[k].L,
                             2 * mA, 2 * pArrB[j].L, 2 * pArrC[k].L);

            for (int iA = 0; iA < canoBasisA[pArrA[i].L - mA_abs].size(); iA++)
            {
              for (int iB = 0; iB < canoBasisB[j].size(); iB++)
              {
                for (int iC = 0; iC < canoBasisC[k].size(); iC++)
                {
                  int n12A = gArrA[canoBasisA[pArrA[i].L - mA_abs][iA]].n12(),
                      n22A = gArrA[canoBasisA[pArrA[i].L - mA_abs][iA]].n22(),
                      n12B = gArrB[canoBasisB[j][iB]].n12(),
                      n22B = gArrB[canoBasisB[j][iB]].n22(),
                      n12C = gArrC[canoBasisC[k][iC]].n12(),
                      n22C = gArrC[canoBasisC[k][iC]].n22();

                  //auto acc =
                  //    canoCoeffs.access(n12A, n22A, n12B, n22B, n12C, n22C);

                  auto [allow, canoCoeffs] = wigU3U2db.get(
                    pArrA[0].toprow.n13(), pArrA[0].toprow.n23(), 
                    pArrA[0].toprow.n33(), n12A, n22A,
                    pArrB[0].toprow.n13(), pArrB[0].toprow.n23(), 
                    pArrB[0].toprow.n33(), n12B, n22B,
                    pArrC[0].toprow.n13(), pArrC[0].toprow.n23(), 
                    pArrC[0].toprow.n33(), n12C, n22C);

                  int jj1 =
                          gArrA[canoBasisA[pArrA[i].L - mA_abs][iA]].Lambdax2(),
                      mm1 = gArrA[canoBasisA[pArrA[i].L - mA_abs][iA]]
                                .MLambdax2(),
                      jj2 = gArrB[canoBasisB[j][iB]].Lambdax2(),
                      mm2 = gArrB[canoBasisB[j][iB]].MLambdax2(),
                      jj3 = gArrC[canoBasisC[k][iC]].Lambdax2(),
                      mm3 = gArrC[canoBasisC[k][iC]].MLambdax2();
                  if (mA < 0)
                    mm1 *= -1;

                  if ((allow == false) || ((mm1 + mm2) != mm3))
                    continue;

                  auto wigU3U2U1 = canoCoeffs *
                                cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

                  int phase =
                      ((pArrA[i].toprow.lam() + pArrA[i].L +
                        gArrA[canoBasisA[pArrA[i].L - mA_abs][iA]].n22()) %
                           2 == 0) ? 1 : -1;

                  for (int r = 0; r < rhomax; r++)
                    for (int k1 = 0; k1 < pArrA[i].kappa; k1++)
                      for (int k2 = 0; k2 < pArrB[j].kappa; k2++)
                        for (int k3 = 0; k3 < pArrC[k].kappa; k3++)
                          if (mA >= 0)
                            coeffs_temp((k1 * pArrB[j].kappa + k2) *
                                         pArrC[k].kappa + k3, r) +=
                                transformationA[pArrA[i].L - mA_abs](iA, k1) *
                                transformationB[j](iB, k2) *
                                transformationC[k](iC, k3) * wigU3U2U1(r);
                          else
                            coeffs_temp((k1 * pArrB[j].kappa + k2) *
                                         pArrC[k].kappa + k3, r) += phase *
                                transformationA[pArrA[i].L - mA_abs](iA, k1) *
                                transformationB[j](iB, k2) *
                                transformationC[k](iC, k3) * wigU3U2U1(r);
                }
              }
            }
            coeffs_temp /= cgc1;
            coeffs.push_back(coeffs_temp);
          }
        }
      }
    }

    return coeffs;
  }

#define U3LIB_ENABLE_WIGU3SO3_CACHE
  class WigU3SO3_db
  {
#ifdef U3LIB_ENABLE_WIGU3SO3_CACHE
    using key_type = std::array<int, 9>;
    std::map<key_type, Eigen::MatrixXd> map_;
#endif

    qcx_utils::SO3_CGC& cgcs_;
    u3::WigU3U2_db& wigU3U2db_;

    auto calculate(const irrep& u3a, int la, 
                   const irrep& u3b, int lb, 
                   const irrep& u3c, int lc, 
                   int rhomax)
    {
      std::vector<u3::physicalvector> 
        pArrA(1, u3::physicalvector(u3a, la)),
        pArrB(1, u3::physicalvector(u3b, lb)), 
        pArrC(1, u3::physicalvector(u3c, lc)); 

      std::vector<std::tuple<int, int, int>> tArr;

      /* return one coefficient matrix
        row = (ka * max(kb) + kb) * max(kc) + kc
        column = rho
      */
      return WigU3SO3(rhomax, pArrA, pArrB, pArrC, tArr, cgcs_, wigU3U2db_)[0];
    }

  public:
    WigU3SO3_db(qcx_utils::SO3_CGC& cgcs, u3::WigU3U2_db& wigU3U2db)
        : cgcs_(cgcs), wigU3U2db_(wigU3U2db)
    { }

#ifdef U3LIB_ENABLE_WIGU3SO3_CACHE
    WigU3SO3_db(const WigU3SO3_db&) = delete;
    WigU3SO3_db& operator=(const WigU3SO3_db&) = delete;
#endif

    auto get(int lam1, int mu1, int l1, 
             int lam2, int mu2, int l2,
             int lam3, int mu3, int l3)
    {
      // normalized irreps:
      irrep u3a(lam1, mu1);
      irrep u3b(lam2, mu2);
      irrep u3c(lam3, mu3);
      int rhomax = Rho(u3a, u3b, u3c);
      int kappamaxa = Kappa(u3a,l1);
      int kappamaxb = Kappa(u3b,l2);
      int kappamaxc = Kappa(u3c,l3);

      bool flag = false;
      Eigen::MatrixXd m;

      if (rhomax > 0 && kappamaxa > 0 && kappamaxb > 0 && kappamaxc > 0 &&
          std::fabs(l1 - l2) <= l3 && (l1 + l2) >= l3)
      {
#ifdef U3LIB_ENABLE_WIGU3SO3_CACHE
        key_type key = {lam1, mu1, l1,
                        lam2, mu2, l2,
                        lam3, mu3, l3};

        auto iter = map_.find(key);
        if (iter == map_.end())
        {
          auto p = map_.try_emplace(
            key, calculate(u3a, l1, u3b, l2, u3c, l3, rhomax));
          assert(p.second);
          iter = p.first;
        }

        m = iter->second;
#else
        m = calculate(u3a, l1, u3b, l2, u3c, l3, rhomax);
#endif
        flag = true;
        
      }

      return std::make_pair(flag, m);
    }
  };
} // namespace u3

#endif
