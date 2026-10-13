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

#include "Definitions.h"
#ifdef TIME_WIGU3SO3
#include <chrono>
#endif

#include "U3basis.h"
#include "U3SO3_to_U3U2U1.h"
#include "U3U2xU3U2.h"

#include <Cache.h>
#include <Int.h>
#include <so3Coeffs/SO3_CGC.h>

namespace u3
{
  // Generate Wigner U3->SO(3) coefficients for multiple values of l1, l2, l3
  inline auto
  WigU3SO3(int rhomax, 
           const std::vector<u3::physicalvector>& pArrA,
           const std::vector<u3::physicalvector>& pArrB,
           const std::vector<u3::physicalvector>& pArrC)
  {
    auto& cgcs = utils::cgcs();
    auto& wigU3U2db = u3::wigU3U2db();

    // build Gelfand basis
    std::vector<u3::gelfand> gArrA = u3::gelfandBasis(pArrA[0].irr()),
                             gArrB = u3::gelfandBasis(pArrB[0].irr()),
                             gArrC = u3::gelfandBasis(pArrC[0].irr());

    // the U3xU2 table is keyed by the irrep triple only -- fetch it once for
    // the whole function and resolve rows locally
    auto canoCoeffsTable =
        wigU3U2db.get_coeffs(pArrA[0].irr(), pArrB[0].irr(), pArrC[0].irr());
    assert(canoCoeffsTable);

    // scratch row scaled by the U(1) CGC (rhomax is constant here)
    Eigen::RowVectorXd wigU3U2U1(rhomax);

    // find transformation for physical states belong to the 3rd irrep
    std::vector<std::vector<int>> canoBasisC(pArrC.size());
    std::vector<Eigen::MatrixXd> transformationC(pArrC.size());
    for (int i = 0; i < pArrC.size(); i++)
      std::tie(canoBasisC[i], transformationC[i]) =
          u3::phys_to_cano(pArrC[i], pArrC[i].L(), gArrC);

    // arrays to store the results
    std::vector<U3SO3<coeffs_t>> coeffs;
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

            U3SO3<coeffs_t> coeffs_temp(KA, KB, KC, rhomax);
            
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

                  auto acc = canoCoeffsTable->find_row(gA.n12(), gA.n22(),
                                                       gB.n12(), gB.n22(),
                                                       gC.n12(), gC.n22());

                  int jj1 = gA.Lambdax2(), mm1 = gA.MLambdax2(),
                      jj2 = gB.Lambdax2(), mm2 = gB.MLambdax2(),
                      jj3 = gC.Lambdax2(), mm3 = gC.MLambdax2();

                  if (!acc.allowed || ((mm1 + mm2) != mm3))
                    continue;

                  wigU3U2U1 = canoCoeffsTable->get_row(acc) *
                              cgcs.get({jj1, jj2, jj3, mm1, mm2, mm3});

                  for (int r = 0; r < rhomax; r++)
                    for (int k1 = 0; k1 < KA; k1++)
                      for (int k2 = 0; k2 < KB; k2++)
                        for (int k3 = 0; k3 < KC; k3++)
                          coeffs_temp(k1, k2, k3, r) +=
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

            U3SO3<coeffs_t> coeffs_temp(KA, KB, KC, rhomax);

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

                  // find_row expects normalized labels (n33 subtracted)
                  auto acc = canoCoeffsTable->find_row(
                      n12A - pArrA[0].n33(), n22A - pArrA[0].n33(),
                      n12B - pArrB[0].n33(), n22B - pArrB[0].n33(),
                      n12C - pArrC[0].n33(), n22C - pArrC[0].n33());

                  int jj1 = gArrA[canoBasisA[LA - mA_abs][iA]].Lambdax2(),
                      mm1 = gArrA[canoBasisA[LA - mA_abs][iA]].MLambdax2(),
                      jj2 = gArrB[canoBasisB[j][iB]].Lambdax2(),
                      mm2 = gArrB[canoBasisB[j][iB]].MLambdax2(),
                      jj3 = gArrC[canoBasisC[k][iC]].Lambdax2(),
                      mm3 = gArrC[canoBasisC[k][iC]].MLambdax2();
                  if (mA < 0) mm1 *= -1;

                  if (!acc.allowed || ((mm1 + mm2) != mm3)) continue;

                  wigU3U2U1 = canoCoeffsTable->get_row(acc) *
                              cgcs.get({jj1, jj2, jj3, mm1, mm2, mm3});

                  int phase =
                      ((pArrA[i].lam() + LA +
                        gArrA[canoBasisA[LA - mA_abs][iA]].n22()) %
                           2 == 0) ? 1 : -1;

                  for (int r = 0; r < rhomax; r++)
                    for (int k1 = 0; k1 < KA; k1++)
                      for (int k2 = 0; k2 < KB; k2++)
                        for (int k3 = 0; k3 < KC; k3++)
                          if (mA >= 0)
                            coeffs_temp(k1, k2, k3, r) +=
                                transformationA[LA - mA_abs](iA, k1) *
                                transformationB[j](iB, k2) *
                                transformationC[k](iC, k3) * wigU3U2U1(r);
                          else
                            coeffs_temp(k1, k2, k3, r) += phase *
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

  // Database of Wigner U3->SO(3) coefficients
  class WigU3SO3_db
  {
    // using key_type = std::array<uint8_t, 9>;
    using key_type = std::pair<uint64_t, uint8_t>;

    using value_type = coeffs_t;

#ifdef U3LIB_ENABLE_WIGU3SO3_CACHE
    utils::Cache<key_type, value_type> cache_;
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
#ifdef TIME_WIGU3SO3
      static std::chrono::duration<double> wigu3so3_time(0);
      static double wigu3so3_last_print = 0.0;
      auto start = std::chrono::high_resolution_clock::now();
#endif

      auto& cgcs = utils::cgcs();
      auto& wigU3U2db = u3::wigU3U2db();

      // build Gelfand basis
      std::vector<u3::gelfand> gArrA = u3::gelfandBasis(pA.irr()),
                               gArrB = u3::gelfandBasis(pB.irr()),
                               gArrC = u3::gelfandBasis(pC.irr());

      int LA = pA.L(), LB = pB.L(), LC = pC.L(), KA = pA.kappa(),
          KB = pB.kappa(), KC = pC.kappa();

      // The U3xU2 cache key is the irrep triple, constant for the whole
      // Gelfand loop nest below: fetch the coefficient table once and resolve
      // rows locally (no per-iteration concurrent-map visits).
      auto canoCoeffs =
          wigU3U2db.get_coeffs(pA.irr(), pB.irr(), pC.irr());
      assert(canoCoeffs);   // rhomax > 0 guaranteed by the caller

      // find transformation for the 3rd state
      auto [canoC, transC] = u3::phys_to_cano(pC, LC, gArrC);

      U3SO3<coeffs_t> coeffs(KA, KB, KC, rhomax);

      // scratch row scaled by the U(1) CGC; preallocated once (same rhomax
      // for every row), preserving the exact multiplication order of the
      // previous per-iteration `canoCoeffs * cgc` temporary
      Eigen::RowVectorXd wigU3U2U1(rhomax);

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

            int jj1 = gA.Lambdax2(), mm1 = gA.MLambdax2(),
                jj2 = gB.Lambdax2(), mm2 = gB.MLambdax2(),
                jj3 = gC.Lambdax2(), mm3 = gC.MLambdax2();

            auto cgc = cgcs.get({jj1, jj2, jj3, mm1, mm2, mm3});
            if (std::fabs(cgc) == 0.0) continue;

            auto acc = canoCoeffs->find_row(gA.n12(), gA.n22(),
                                            gB.n12(), gB.n22(),
                                            gC.n12(), gC.n22());
            if (!acc.allowed) continue;

            wigU3U2U1 = canoCoeffs->get_row(acc) * cgc;

            for (int k1 = 0; k1 < KA; k1++)
            {
              double a = transA(iA, k1);
              for (int k2 = 0; k2 < KB; k2++)
              {
                double ab = a * transB(iB, k2);
                for (int k3 = 0; k3 < KC; k3++)
                {
                  double abc = ab * transC(iC, k3);

                  for (int r = 0; r < rhomax; r++)
                    coeffs(k1, k2, k3, r) += abc * wigU3U2U1(r);
                }
              }
            }
          }
        }
      }

      for (size_t i = 0; i < coeffs.size(); i++)
      {
        if (std::fabs(coeffs[i]) < 1e-12) coeffs[i] = 0.0;
        else coeffs[i] /= cgc1;
      }

#ifdef TIME_WIGU3SO3
      auto end = std::chrono::high_resolution_clock::now();
      wigu3so3_time += end - start; // duration in seconds

      // Print occasionally
      if (wigu3so3_time.count() - wigu3so3_last_print >= 10.0)
      {
        std::cout << "Accumulated WigU3SO3 time: " << wigu3so3_time.count()
                  << " s\n";
        wigu3so3_last_print = wigu3so3_time.count();
      }
#endif

      return coeffs.data();
    }

  public:
    U3SO3<value_type> get_U3SO3(int lam1, int mu1, int l1, int lam2,
                                int mu2, int l2, int lam3, int mu3, int l3)
    {
      // normalized irreps:
      u3::irrep u3a(lam1, mu1), u3b(lam2, mu2), u3c(lam3, mu3);
      int rhomax = u3::Rho(u3a, u3b, u3c);
      u3::physicalvector pA(u3a, l1), pB(u3b, l2), pC(u3c, l3);

      if ((rhomax == 0) || (pA.kappa() == 0) || (pB.kappa() == 0) ||
          (pC.kappa() == 0) || (u2::mult(l1, l2, l3) == 0))
        return {};

#ifdef U3LIB_ENABLE_WIGU3SO3_CACHE
      std::array<utils::UInt8, 9> safe_key = {lam1, mu1,  l1,  lam2, mu2,
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

      return {std::move(*p_result_), 
              static_cast<size_t>(pA.kappa()),
              static_cast<size_t>(pB.kappa()),
              static_cast<size_t>(pC.kappa()), 
              static_cast<size_t>(rhomax)};
    }

    value_type get(int lam1, int mu1, int l1, int lam2,
                   int mu2, int l2, int lam3, int mu3, int l3)
    {
      // normalized irreps:
      u3::irrep u3a(lam1, mu1), u3b(lam2, mu2), u3c(lam3, mu3);
      int rhomax = u3::Rho(u3a, u3b, u3c);
      u3::physicalvector pA(u3a, l1), pB(u3b, l2), pC(u3c, l3);

      if ((rhomax == 0) || (pA.kappa() == 0) || (pB.kappa() == 0) ||
          (pC.kappa() == 0) || (u2::mult(l1, l2, l3) == 0))
        return {};

#ifdef U3LIB_ENABLE_WIGU3SO3_CACHE
      std::array<utils::UInt8, 9> safe_key = {lam1, mu1,  l1,  lam2, mu2,
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

      return *p_result_;
    }

    friend WigU3SO3_db& wigU3SO3db();

#ifdef U3LIB_ENABLE_WIGU3SO3_CACHE
    void clear() { cache_.clear(); }
#endif
  };

  inline WigU3SO3_db& wigU3SO3db()
  {
    static WigU3SO3_db wigU3SO3db;

    return wigU3SO3db;
  }

} // namespace u3

#endif
