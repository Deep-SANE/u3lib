#ifndef U3LIB_U3SO3xU3SO3_H
#define U3LIB_U3SO3xU3SO3_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "U3basis.h"
#include "U3SO3_to_U3U2U1.h"
#include "U3U2xU3U2.h"

#include <qcx_Cache.h>
#include <qcx_Int.h>
#include <so3Coeffs/SO3_CGC.h>

namespace u3
{
  // generate Wigner U3->SO(3) coefficients for multiple values of l1, l2, l3
  inline auto
  WigU3SO3(int rhomax, 
           const std::vector<u3::physicalvector>& pArrA,
           const std::vector<u3::physicalvector>& pArrB,
           const std::vector<u3::physicalvector>& pArrC)
  {
    auto& cgcs = qcx_utils::cgcs();
    auto& wigU3U2db = u3::wigU3U2db();

    // build Gelfand basis
    std::vector<u3::gelfand> gArrA = u3::gelfandBasis(pArrA[0].irr()),
                             gArrB = u3::gelfandBasis(pArrB[0].irr()),
                             gArrC = u3::gelfandBasis(pArrC[0].irr());

    // find transformation for physical states belong to the 3rd irrep
    std::vector<std::vector<int>> canoBasisC(pArrC.size());
    std::vector<Eigen::MatrixXd> transformationC(pArrC.size());
    for (int i = 0; i < pArrC.size(); i++)
      std::tie(canoBasisC[i], transformationC[i]) =
          u3::phys_to_cano(pArrC[i], pArrC[i].L(), gArrC);

    std::vector<Eigen::MatrixXd> coeffs;
    std::vector<std::tuple<int, int, int>> tArr;

#ifdef ADD_SU3LIB_PHASE_U3SO3
    // find transformation for physical states belong to the 1st irrep
    std::vector<std::vector<int>> canoBasisA(pArrA.size());
    std::vector<Eigen::MatrixXd> transformationA(pArrA.size());
    for (int i = 0; i < pArrA.size(); i++)
      std::tie(canoBasisA[i], transformationA[i]) =
          u3::phys_to_cano(pArrA[i], pArrA[i].L(), gArrA);

    int i, j, k;
    for (j = 0; j < pArrB.size(); j++)
    {
      auto LB = pArrB[j].L();
      auto KB = pArrB[j].kappa();
      auto [canoBasisB, transformationB] = u3::phys_to_cano_allM(pArrB[j], gArrB);

      for (i = 0; i < pArrA.size(); i++)
      {
        auto LA = pArrA[i].L();
        auto KA = pArrA[i].kappa();
        const auto& canoA = canoBasisA[i];
        const auto& transA = transformationA[i];

        for (k = 0; k < pArrC.size(); k++)
        {
          auto LC = pArrC[k].L();
          auto KC = pArrC[k].kappa();
          const auto& canoC = canoBasisC[k];
          const auto& transC = transformationC[k];

          if (std::fabs(LA - LB) <= LC && (LA + LB) >= LC)
          {
            tArr.emplace_back(i, j, k);

            Eigen::MatrixXd coeffs_temp(KA * KB * KC, rhomax);
            coeffs_temp.setZero();

            int mB = LC - LA;
            const auto& canoB = canoBasisB[LB - mB];
            const auto& transB = transformationB[LB - mB];

            double cgc1 =
                cgcs.get({2 * LA, 2 * LB, 2 * LC, 2 * LA, 2 * mB, 2 * LC});

            for (int iA = 0; iA < canoA.size(); iA++)
            {
              const auto& gA = gArrA[canoA[iA]];

              for (int iB = 0; iB < canoB.size(); iB++)
              {
                const auto& gB = gArrB[canoB[iB]];

                for (int iC = 0; iC < canoC.size(); iC++)
                {
                  const auto& gC = gArrC[canoC[iC]];

                  auto [allow, canoCoeffs] = wigU3U2db.get(
                    gA.n13(), gA.n23(), gA.n33(), gA.n12(), gA.n22(),
                    gB.n13(), gB.n23(), gB.n33(), gB.n12(), gB.n22(),
                    gC.n13(), gC.n23(), gC.n33(), gC.n12(), gC.n22());

                  int jj1 = gA.Lambdax2(), mm1 = gA.MLambdax2(),
                      jj2 = gB.Lambdax2(), mm2 = gB.MLambdax2(),
                      jj3 = gC.Lambdax2(), mm3 = gC.MLambdax2();

                  if ((allow == false) || ((mm1 + mm2) != mm3))
                    continue;

                  auto wigU3U2U1 =
                      canoCoeffs * cgcs.get({jj1, jj2, jj3, mm1, mm2, mm3});

                  for (int r = 0; r < rhomax; r++)
                    for (int k1 = 0; k1 < KA; k1++)
                      for (int k2 = 0; k2 < KB; k2++)
                        for (int k3 = 0; k3 < KC; k3++)
                          coeffs_temp((k1 * KB + k2) * KC + k3, r) +=
                                       transA(iA, k1) * transB(iB, k2) *
                                       transC(iC, k3) * wigU3U2U1(r);
                }
              }
            }
            coeffs_temp /= cgc1;
            coeffs.push_back(coeffs_temp);
          }
        }
      }
    }
#else
    // find transformation for physical states belong to the 2nd irrep
    std::vector<std::vector<int>> canoBasisB(pArrB.size());
    std::vector<Eigen::MatrixXd> transformationB(pArrB.size());
    for (int i = 0; i < pArrB.size(); i++)
      std::tie(canoBasisB[i], transformationB[i]) =
          u3::phys_to_cano(pArrB[i], pArrB[i].L(), gArrB);

    int i, j, k;
    for (i = 0; i < pArrA.size(); i++)
    {
      auto LA = pArrA[i].L();
      auto KA = pArrA[i].kappa();

      std::vector<std::vector<int>> canoBasisA(1);
      std::vector<Eigen::MatrixXd> transformationA(1);
        std::tie(canoBasisA[0], transformationA[0]) =
          u3::phys_to_cano(pArrA[i], LA, gArrA);

      for (j = 0; j < pArrB.size(); j++)
      {
        auto LB = pArrB[j].L();
        auto KB = pArrB[j].kappa();

        for (k = 0; k < pArrC.size(); k++)
        {
          auto LC = pArrC[k].L();
          auto KC = pArrC[k].kappa();

          if (std::fabs(LA - LB) <= LC && (LA + LB) >= LC)
          {
            tArr.emplace_back(i, j, k);

            Eigen::MatrixXd coeffs_temp(KA * KB * KC, rhomax);
            coeffs_temp.setZero();

            int mA = LC - LB, mA_abs;
            if (mA >= 0) mA_abs = mA;
            else mA_abs = -mA;

            while ((LA - mA_abs + 1) > canoBasisA.size())
            {
              int m_new = LA - canoBasisA.size();
              std::vector<int> canoBasisA_new =
                  u3::expand_to_gelfand(m_new, gArrA);

              Eigen::MatrixXd transformationA_new = u3::lowerM(KA, 
                LA, m_new, gArrA, transformationA.back(), 
                canoBasisA.back(), canoBasisA_new);

              canoBasisA.push_back(canoBasisA_new);
              transformationA.push_back(transformationA_new);
            }

            double cgc1 =
                cgcs.get({2 * LA, 2 * LB, 2 * LC, 2 * mA, 2 * LB, 2 * LC});

            for (int iA = 0; iA < canoBasisA[LA - mA_abs].size(); iA++)
            {
              for (int iB = 0; iB < canoBasisB[j].size(); iB++)
              {
                for (int iC = 0; iC < canoBasisC[k].size(); iC++)
                {
                  int n12A = gArrA[canoBasisA[LA - mA_abs][iA]].n12(),
                      n22A = gArrA[canoBasisA[LA - mA_abs][iA]].n22(),
                      n12B = gArrB[canoBasisB[j][iB]].n12(),
                      n22B = gArrB[canoBasisB[j][iB]].n22(),
                      n12C = gArrC[canoBasisC[k][iC]].n12(),
                      n22C = gArrC[canoBasisC[k][iC]].n22();

                  //auto acc =
                  //    canoCoeffs.access(n12A, n22A, n12B, n22B, n12C, n22C);

                  auto [allow, canoCoeffs] = wigU3U2db.get(
                    pArrA[0].n13(), pArrA[0].n23(), pArrA[0].n33(), n12A, n22A,
                    pArrB[0].n13(), pArrB[0].n23(), pArrB[0].n33(), n12B, n22B,
                    pArrC[0].n13(), pArrC[0].n23(), pArrC[0].n33(), n12C, n22C);

                  int jj1 = gArrA[canoBasisA[LA - mA_abs][iA]].Lambdax2(),
                      mm1 = gArrA[canoBasisA[LA - mA_abs][iA]].MLambdax2(),
                      jj2 = gArrB[canoBasisB[j][iB]].Lambdax2(),
                      mm2 = gArrB[canoBasisB[j][iB]].MLambdax2(),
                      jj3 = gArrC[canoBasisC[k][iC]].Lambdax2(),
                      mm3 = gArrC[canoBasisC[k][iC]].MLambdax2();
                  if (mA < 0) mm1 *= -1;

                  if ((allow == false) || ((mm1 + mm2) != mm3)) continue;

                  auto wigU3U2U1 =
                      canoCoeffs * cgcs.get({jj1, jj2, jj3, mm1, mm2, mm3});

                  int phase =
                      ((pArrA[i].lam() + LA +
                        gArrA[canoBasisA[LA - mA_abs][iA]].n22()) %
                           2 == 0) ? 1 : -1;

                  for (int r = 0; r < rhomax; r++)
                    for (int k1 = 0; k1 < KA; k1++)
                      for (int k2 = 0; k2 < KB; k2++)
                        for (int k3 = 0; k3 < KC; k3++)
                          if (mA >= 0)
                            coeffs_temp((k1 * KB + k2) * KC + k3, r) +=
                                transformationA[LA - mA_abs](iA, k1) *
                                transformationB[j](iB, k2) *
                                transformationC[k](iC, k3) * wigU3U2U1(r);
                          else
                            coeffs_temp((k1 * KB + k2) * KC + k3, r) += phase *
                                transformationA[LA - mA_abs](iA, k1) *
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
#endif

    return std::make_pair(tArr, coeffs);
  }

#define U3LIB_ENABLE_WIGU3SO3_CACHE

  class WigU3SO3_db
  {
    // using key_type = std::array<uint8_t, 9>;
    using key_type = std::pair<uint64_t, uint8_t>;
    using value_type = Eigen::MatrixXd;

#ifdef U3LIB_ENABLE_WIGU3SO3_CACHE
    qcx_utils::Cache<key_type, value_type> cache_;
#endif

    inline static value_type* p_result_ = nullptr;
#pragma omp threadprivate(p_result_)

  public:
    static void initialize()
    {
      assert(p_result_ == nullptr);
      p_result_ = new value_type;
    }

    static void finalize()
    {
      assert(p_result_ != nullptr);
      delete p_result_;
      p_result_ = nullptr;
    }

    static void initialize_thread()
    {
#ifdef USE_OPENMP
      if (omp_get_thread_num() == 0)
      {
        assert(p_result_ != nullptr);
      }
      else
      {
        assert(p_result_ == nullptr);
        p_result_ = new value_type;
      }
#endif
    }

    static void finalize_thread()
    {
#ifdef USE_OPENMP
      if (omp_get_thread_num() != 0)
      {
        assert(p_result_ != nullptr);
        delete p_result_;
        p_result_ = nullptr;
      }
#endif
    }

  private:
    // This class is a singleton, and cannot be instantiated directly:
    WigU3SO3_db() = default;
    WigU3SO3_db(const WigU3SO3_db&) = delete;
    WigU3SO3_db& operator=(const WigU3SO3_db&) = delete;

    // generate Wigner U3->SO(3) coefficients for specific l1, l2, l3
    static value_type compute(int rhomax, const u3::physicalvector& pA,
                              const u3::physicalvector& pB,
                              const u3::physicalvector& pC)
    {
      auto& cgcs = qcx_utils::cgcs();
      auto& wigU3U2db = u3::wigU3U2db();

      // build Gelfand basis
      std::vector<u3::gelfand> gArrA = u3::gelfandBasis(pA.irr()),
                               gArrB = u3::gelfandBasis(pB.irr()),
                               gArrC = u3::gelfandBasis(pC.irr());

      int LA = pA.L(), LB = pB.L(), LC = pC.L(), KA = pA.kappa(),
          KB = pB.kappa(), KC = pC.kappa();

      // find transformation for the 3rd state
      auto [canoC, transC] = u3::phys_to_cano(pC, LC, gArrC);

      Eigen::MatrixXd coeffs(KA * KB * KC, rhomax);
      coeffs.setZero();

#ifdef ADD_SU3LIB_PHASE_U3SO3
      // find transformation for the 1st state
      auto [canoA, transA] = u3::phys_to_cano(pA, LA, gArrA);

      // find transformation for the 2nd state
      int mB = LC - LA;
      auto [canoB, transB] = u3::phys_to_cano(pB, mB, gArrB);

      double cgc1 = cgcs.get({2 * LA, 2 * LB, 2 * LC, 2 * LA, 2 * mB, 2 * LC});
#else
      // find transformation for the 1st state
      int mA = LC - LB;
      auto [canoA, transA] = u3::phys_to_cano(pA, mA, gArrA);

      // find transformation for the 2nd state
      auto [canoB, transB] = u3::phys_to_cano(pB, LB, gArrB);

      double cgc1 = cgcs.get({2 * LA, 2 * LB, 2 * LC, 2 * mA, 2 * LB, 2 * LC});
#endif

      for (int iA = 0; iA < canoA.size(); iA++)
      {
        const auto& gA = gArrA[canoA[iA]];

        for (int iB = 0; iB < canoB.size(); iB++)
        {
          const auto& gB = gArrB[canoB[iB]];

          for (int iC = 0; iC < canoC.size(); iC++)
          {
            const auto& gC = gArrC[canoC[iC]];

            auto [allow, canoCoeffs] = wigU3U2db.get(
              gA.n13(), gA.n23(), gA.n33(), gA.n12(), gA.n22(),
              gB.n13(), gB.n23(), gB.n33(), gB.n12(), gB.n22(),
              gC.n13(), gC.n23(), gC.n33(), gC.n12(), gC.n22());

            int jj1 = gA.Lambdax2(), mm1 = gA.MLambdax2(),
                jj2 = gB.Lambdax2(), mm2 = gB.MLambdax2(),
                jj3 = gC.Lambdax2(), mm3 = gC.MLambdax2();

            if ((allow == false) || ((mm1 + mm2) != mm3)) continue;

            auto wigU3U2U1 =
                canoCoeffs * cgcs.get({jj1, jj2, jj3, mm1, mm2, mm3});

            for (int r = 0; r < rhomax; r++)
              for (int k1 = 0; k1 < KA; k1++)
                for (int k2 = 0; k2 < KB; k2++)
                  for (int k3 = 0; k3 < KC; k3++)
                    coeffs((k1 * KB + k2) * KC + k3, r) 
                      += transA(iA, k1) * transB(iB, k2) *
                         transC(iC, k3) * wigU3U2U1(r);
          }
        }
      }

      coeffs /= cgc1;

      for (size_t i = 0; i < coeffs.rows(); i++)
        for (size_t j = 0; j < coeffs.cols(); j++)
          if (std::fabs(coeffs(i, j)) < 1e-12)
            coeffs(i, j) = 0.0;
    
      return coeffs;
    }

  public:
    std::pair<bool, value_type> get(int lam1, int mu1, int l1, int lam2,
                                    int mu2, int l2, int lam3, int mu3, int l3)
    {
      // normalized irreps:
      u3::irrep u3a(lam1, mu1), u3b(lam2, mu2), u3c(lam3, mu3);
      int rhomax = u3::Rho(u3a, u3b, u3c);
      u3::physicalvector pA(u3a, l1), pB(u3b, l2), pC(u3c, l3);

      if ((rhomax == 0) || (pA.kappa() == 0) || (pB.kappa() == 0) ||
          (pC.kappa() == 0) || (u2::mult(l1, l2, l3) == 0))
        return {false, value_type()};

#ifdef U3LIB_ENABLE_WIGU3SO3_CACHE
      std::array<qcx_utils::UInt8, 9> safe_key = {lam1, mu1,  l1,  lam2, mu2,
                                                  l2,   lam3, mu3, l3};

      uint64_t key_1 = 
        (uint64_t(safe_key[0]) << 56) | (uint64_t(safe_key[1]) << 48) |
        (uint64_t(safe_key[2]) << 40) | (uint64_t(safe_key[3]) << 32) |
        (uint64_t(safe_key[4]) << 24) | (uint64_t(safe_key[5]) << 16) |
        (uint64_t(safe_key[6]) << 8) | uint64_t(safe_key[7]);

      key_type key = {key_1, safe_key[8]};
      /*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5],
                      safe_key[6], safe_key[7], safe_key[8]};
      */
      auto producer = [rhomax, &pA, &pB, &pC](key_type)
      { return compute(rhomax, pA, pB, pC); };

      auto accessor = [](const auto& element) { *p_result_ = element.second; };

      cache_.visit(key, producer, accessor);
#else
      *p_result_ = compute(rhomax, pA, pB, pC);
#endif

      return {true, *p_result_};
    }

    friend WigU3SO3_db& wigU3SO3db();
  };

  inline WigU3SO3_db& wigU3SO3db()
  {
    static WigU3SO3_db wigU3SO3db;

    return wigU3SO3db;
  }

} // namespace u3

#endif
