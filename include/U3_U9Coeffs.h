#ifndef U3LIB_U3_U9COEFFS_H
#define U3LIB_U3_U9COEFFS_H 

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <map>
#include <span>
#include <utility>
#include <vector>

#ifdef TIME_U3U9
#include <chrono>
#endif

#include <Eigen/Dense>

#include "U3decomposition.h"
#include "U3_U6Z6Coeffs.h"
#include "U3basis.h"

#include <qcx_ND.h>
#include <so3Coeffs/SO3_CGC.h>
#include <so3Coeffs/SO3_Racah.h>

namespace u3
{
  // input by users: u3a, u3b, u3c, u3d, u3ab, u3cd, u3ac, u3bd, u3abcd
  // all normalized
  inline auto
  U9(const u3::irrep& u3A, const std::vector<u3::gelfand>& gArrA,
     const u3::irrep& u3B, const std::vector<u3::gelfand>& gArrB,
     const u3::irrep& u3C, const std::vector<u3::gelfand>& gArrC,
     const u3::irrep& u3D, const std::vector<u3::gelfand>& gArrD,
     const u3::irrep& u3AB, const std::vector<u3::gelfand>& gArrAB,
     const u3::irrep& u3CD, const std::vector<u3::gelfand>& gArrCD,
     const u3::irrep& u3AC, const std::vector<u3::gelfand>& gArrAC,
     const u3::irrep& u3BD, const std::vector<u3::gelfand>& gArrBD,
     const u3::irrep& u3ABCD, const std::vector<u3::gelfand>& gArrABCD,
     int rhoABmax, int rhoCDmax, int rhoAB_CDmax, int rhoACmax, int rhoBDmax,
     int rhoAC_BDmax)
  {
#ifdef TIME_U3U9
    static std::chrono::duration<double> u3u9_time(0);
    static double u3u9_last_print = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
#endif

    auto& cgcs = qcx_utils::cgcs();
    auto& racahs = qcx_utils::racahs();
    auto& wigU3U2db = u3::wigU3U2db();
    auto& u3u6db = u3::u3u6db();
    auto& u3z6db = u3::u3z6db();
/*
    std::vector<double> U(rhoACmax * rhoBDmax * rhoAB_CDmax * rhoABmax *
                              rhoCDmax * rhoAC_BDmax, 0.0);
*/
    U3U9<coeffs_t> U(rhoACmax, rhoBDmax, rhoAB_CDmax, 
                     rhoABmax, rhoCDmax, rhoAC_BDmax);

    // couple u3AB and u3C
    auto [u3AB_C, rhomaxsAB_C] = u3::decomposition(u3AB, u3C);
    
    for (int i = 0; i < u3AB_C.size(); i++)
    {
      int rhomaxAC_B = u3::Rho(u3AC, u3B, u3AB_C[i]);
      if (rhomaxAC_B == 0) continue;

      int rhoABC_Dmax = u3::Rho(u3AB_C[i], u3D, u3ABCD);
      if (rhoABC_Dmax == 0) continue;

#if defined(U3LIB_ENABLE_U3U6_CACHE) && defined(U3LIB_ENABLE_U3Z6_CACHE)
      const auto U6_1 = u3u6db.get_U3U6(u3AC, u3B, u3ABCD, u3D,
                                        u3AB_C[i], rhomaxAC_B, rhoABC_Dmax, 
                                        u3BD, rhoBDmax, rhoAC_BDmax);
      if (U6_1.empty()) continue;

      const auto U6_2 = u3u6db.get_U3U6(u3AB, u3C, u3ABCD, u3D,
                                        u3AB_C[i], rhomaxsAB_C[i], rhoABC_Dmax, 
                                        u3CD, rhoCDmax, rhoAB_CDmax);
      if (U6_2.empty()) continue;
    
      const auto Z6 = u3z6db.get_U3Z6(u3B, u3A, u3AB_C[i], u3C, 
                                      u3AB, rhoABmax, rhomaxsAB_C[i], 
                                      u3AC, rhoACmax, rhomaxAC_B);
      if (Z6.empty()) continue;
#else
      u3::gelfands gArrABC = u3::gelfandBasis(u3AB_C[i]);

      const auto U6_1 =
          u3::U6(u3AC, gArrAC, u3B, gArrB, u3ABCD, gArrABCD, u3D, gArrD,
                 u3AB_C[i], gArrABC, rhomaxAC_B, rhoABC_Dmax, u3BD, gArrBD,
                 rhoBDmax, rhoAC_BDmax);

      const auto U6_2 =
          u3::U6(u3AB, gArrAB, u3C, gArrC, u3ABCD, gArrABCD, u3D, gArrD,
                 u3AB_C[i], gArrABC, rhomaxsAB_C[i], rhoABC_Dmax, u3CD, gArrCD,
                 rhoCDmax, rhoAB_CDmax);
      const auto Z6 =
          u3::Z6(u3B, gArrB, u3A, gArrA, u3AB_C[i], gArrABC, u3C, gArrC, u3AB,
                 gArrAB, rhoABmax, rhomaxsAB_C[i], u3AC, gArrAC, rhoACmax,
                 rhomaxAC_B);
#endif
      
      for (int rAC = 0; rAC < rhoACmax; rAC++)
        for (int rBD = 0; rBD < rhoBDmax; rBD++)
          for (int rAB_CD = 0; rAB_CD < rhoAB_CDmax; rAB_CD++)
            for (int rAB = 0; rAB < rhoABmax; rAB++)
              for (int rCD = 0; rCD < rhoCDmax; rCD++)
                for (int rAC_BD = 0; rAC_BD < rhoAC_BDmax; rAC_BD++)
                  for (int rAC_B = 0; rAC_B < rhomaxAC_B; rAC_B++)
                    for (int rAB_C = 0; rAB_C < rhomaxsAB_C[i]; rAB_C++)
                      for (int rABC_D = 0; rABC_D < rhoABC_Dmax; rABC_D++)
                      {
#if defined(U3LIB_ENABLE_U3U6_CACHE) && defined(U3LIB_ENABLE_U3Z6_CACHE)
                        U(rAC, rBD, rAB_CD, rAB, rCD, rAC_BD)
                        += U6_1(rAC_B, rABC_D, rBD, rAC_BD) 
                         * Z6(rAB, rAB_C, rAC, rAC_B) 
                         * U6_2(rAB_C, rABC_D, rCD, rAB_CD);
#else
                        U(rAC, rBD, rAB_CD, rAB, rCD, rAC_BD) 
                        += U6_1[((rAC_B * rhoABC_Dmax + rABC_D) * rhoBDmax +
                                  rBD) * rhoAC_BDmax + rAC_BD] 
                         * Z6[((rAB * rhomaxsAB_C[i] + rAB_C) * rhoACmax +
                                rAC) * rhomaxAC_B + rAC_B] 
                         * U6_2[((rAB_C * rhoABC_Dmax + rABC_D) * rhoCDmax +
                                  rCD) * rhoAB_CDmax + rAB_CD];
#endif
                      }
    }

    for (int i = 0; i < U.size(); i++)
      if (std::fabs(U[i]) < 1e-14) U[i] = 0.0;

#ifdef TIME_U3U9
    auto end = std::chrono::high_resolution_clock::now();
    u3u9_time += end - start;  // duration in seconds

    // Print occasionally 
    if (u3u9_time.count() - u3u9_last_print >= 1.0) 
    {
        std::cout << "Accumulated U3U9 time: " << u3u9_time.count() << " s\n";
        u3u9_last_print = u3u9_time.count();
    }
#endif

    return U;
  }

  // Database of U9 coefficients
  class U3U9_db
  {
    // using key_type = std::array<uint8_t, 18>;
    using key_type = std::tuple<uint64_t, uint64_t, uint16_t>;
    using value_type = coeffs_t;

#ifdef U3LIB_ENABLE_U3U9_CACHE
    qcx_utils::Cache<key_type, value_type> cache_;
    using return_type = std::span<const double>;
#else
    using return_type = std::vector<double>;
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
    U3U9_db() = default;
    U3U9_db(const U3U9_db&) = delete;
    U3U9_db& operator=(const U3U9_db&) = delete;

    static value_type compute(
      const u3::irrep& u3A, const u3::irrep& u3B,
      const u3::irrep& u3C, const u3::irrep& u3D,
      const u3::irrep& u3AB, const u3::irrep& u3CD,
      const u3::irrep& u3AC, const u3::irrep& u3BD,
      const u3::irrep& u3ABCD, 
      int rhoABmax, int rhoCDmax, int rhoAB_CDmax, 
      int rhoACmax, int rhoBDmax, int rhoAC_BDmax)
    {
      std::vector<u3::gelfand> gArrA = u3::gelfandBasis(u3A),
                               gArrB = u3::gelfandBasis(u3B),
                               gArrC = u3::gelfandBasis(u3C),
                               gArrD = u3::gelfandBasis(u3D),
                               gArrAB = u3::gelfandBasis(u3AB),
                               gArrCD = u3::gelfandBasis(u3CD),
                               gArrAC = u3::gelfandBasis(u3AC),
                               gArrBD = u3::gelfandBasis(u3BD),
                               gArrABCD = u3::gelfandBasis(u3ABCD);

      return U9(u3A, gArrA, u3B, gArrB, u3C, gArrC, u3D, gArrD, u3AB, gArrAB,
                u3CD, gArrCD, u3AC, gArrAC, u3BD, gArrBD, u3ABCD, gArrABCD,
                rhoABmax, rhoCDmax, rhoAB_CDmax, rhoACmax, rhoBDmax,
                rhoAC_BDmax).data();
    }

  public:
    // pass by references to avoid copies
    U3U9<return_type> get_U3U9(
      const u3::irrep& u3A, const u3::irrep& u3B, 
      const u3::irrep& u3C, const u3::irrep& u3D,
      const u3::irrep& u3AB, const u3::irrep& u3CD,
      const u3::irrep& u3AC, const u3::irrep& u3BD, 
      const u3::irrep& u3ABCD,
      int rhoACmax, int rhoBDmax, int rhoABCDmax,
      int rhoABmax, int rhoCDmax, int rhoACBDmax)
    {
      auto coeffs =
          get(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD, rhoABmax,
              rhoCDmax, rhoABCDmax, rhoACmax, rhoBDmax, rhoACBDmax);

      return {std::move(coeffs),
              static_cast<size_t>(rhoACmax),
              static_cast<size_t>(rhoBDmax),
              static_cast<size_t>(rhoABCDmax),
              static_cast<size_t>(rhoABmax),
              static_cast<size_t>(rhoCDmax),
              static_cast<size_t>(rhoACBDmax)};
    }

    // pass by values to normalize irreps
    return_type get(u3::irrep u3A, u3::irrep u3B, u3::irrep u3C, u3::irrep u3D,
                    u3::irrep u3AB, u3::irrep u3CD, u3::irrep u3AC,
                    u3::irrep u3BD, u3::irrep u3ABCD, int rhoABmax,
                    int rhoCDmax, int rhoABCDmax, int rhoACmax, int rhoBDmax,
                    int rhoACBDmax)
    {
      if ((rhoACmax == 0) || (rhoBDmax == 0) || (rhoABCDmax == 0) ||
          (rhoABmax == 0) || (rhoCDmax == 0) || (rhoACBDmax == 0))
        return {};

      u3A.normalize();
      u3B.normalize();
      u3C.normalize();
      u3D.normalize();
      u3AB.normalize();
      u3CD.normalize();
      u3AC.normalize();
      u3BD.normalize();
      u3ABCD.normalize();

      assert(rhoABmax == Rho(u3A, u3B, u3AB));
      assert(rhoCDmax == Rho(u3C, u3D, u3CD));
      assert(rhoABCDmax == Rho(u3AB, u3CD, u3ABCD));
      assert(rhoACmax == Rho(u3A, u3C, u3AC));
      assert(rhoBDmax == Rho(u3B, u3D, u3BD));
      assert(rhoACBDmax == Rho(u3AC, u3BD, u3ABCD));

#ifdef U3LIB_ENABLE_U3U9_CACHE
      std::array<qcx_utils::UInt8, 18> safe_key = {
        u3A.n13(), u3A.n23(), u3B.n13(), u3B.n23(),
        u3C.n13(), u3C.n23(), u3D.n13(), u3D.n23(),
        u3AB.n13(), u3AB.n23(), u3CD.n13(), u3CD.n23(),
        u3AC.n13(), u3AC.n23(), u3BD.n13(), u3BD.n23(),
        u3ABCD.n13(), u3ABCD.n23()};

      uint64_t key_1 = 
        (uint64_t(safe_key[0]) << 56) | (uint64_t(safe_key[1]) << 48) |
        (uint64_t(safe_key[2]) << 40) | (uint64_t(safe_key[3]) << 32) |
        (uint64_t(safe_key[4]) << 24) | (uint64_t(safe_key[5]) << 16) |
        (uint64_t(safe_key[6]) << 8) | uint64_t(safe_key[7]);

      uint64_t key_2 = 
        (uint64_t(safe_key[8]) << 56) | (uint64_t(safe_key[9]) << 48) |
        (uint64_t(safe_key[10]) << 40) | (uint64_t(safe_key[11]) << 32) |
        (uint64_t(safe_key[12]) << 24) | (uint64_t(safe_key[13]) << 16) |
        (uint64_t(safe_key[14]) << 8) | uint64_t(safe_key[15]);

      uint16_t key_3 = 
        (uint16_t(safe_key[16]) << 8) | uint16_t(safe_key[17]);

      key_type key = {key_1, key_2, key_3};
      /*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5],
                      safe_key[6], safe_key[7], safe_key[8],
                      safe_key[9], safe_key[10], safe_key[11],
                      safe_key[12], safe_key[13], safe_key[14],
                      safe_key[15], safe_key[16], safe_key[17]};
      */
      auto producer = [&u3A, &u3B, &u3C, &u3D, &u3AB, &u3CD, &u3AC, &u3BD, &u3ABCD, 
                       &rhoABmax, &rhoCDmax, &rhoABCDmax,
                       &rhoACmax, &rhoBDmax, &rhoACBDmax](key_type)
      { 
        return compute(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD, 
                       rhoABmax, rhoCDmax, rhoABCDmax,
                       rhoACmax, rhoBDmax, rhoACBDmax);
      };

      auto accessor = [](const auto& element) { *p_result_ = element.second; };

      cache_.visit(key, producer, accessor);
#else
      *p_result_ = compute(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD, 
                           rhoABmax, rhoCDmax, rhoABCDmax,
                           rhoACmax, rhoBDmax, rhoACBDmax);
#endif

      return *p_result_;
    }

    // TO BE REMOVED IN THE FUTURE:
    value_type get(
      int n13a, int n23a, int n33a, int n13b, int n23b, int n33b, 
      int n13c, int n23c, int n33c, int n13d, int n23d, int n33d,
      int n13ab, int n23ab, int n33ab, int n13cd, int n23cd, int n33cd,
      int n13ac, int n23ac, int n33ac, int n13bd, int n23bd, int n33bd,
      int n13abcd, int n23abcd, int n33abcd)
    {
      u3::irrep u3A(n13a - n33a, n23a - n33a, 0), 
                u3B(n13b - n33b, n23b - n33b, 0), 
                u3C(n13c - n33c, n23c - n33c, 0), 
                u3D(n13d - n33d, n23d - n33d, 0), 
                u3AB(n13ab - n33ab, n23ab - n33ab, 0), 
                u3CD(n13cd - n33cd, n23cd - n33cd, 0), 
                u3AC(n13ac - n33ac, n23ac - n33ac, 0),
                u3BD(n13bd - n33bd, n23bd - n33bd, 0),
                u3ABCD(n13abcd - n33abcd, n23abcd - n33abcd, 0);

      int rhoABmax = Rho(u3A,u3B,u3AB), rhoCDmax = Rho(u3C,u3D,u3CD),
          rhoACmax = Rho(u3A,u3C,u3AC), rhoBDmax = Rho(u3B,u3D,u3BD),
          rhoABCDmax = Rho(u3AB,u3CD,u3ABCD), 
          rhoACBDmax = Rho(u3AC,u3BD,u3ABCD);

      if (rhoABmax == 0 || rhoCDmax == 0 || rhoACmax == 0 || 
          rhoBDmax == 0 || rhoABCDmax == 0 || rhoACBDmax == 0)
        return {};

#ifdef U3LIB_ENABLE_U3U9_CACHE
      std::array<qcx_utils::UInt8, 18> safe_key = {
        u3A.n13(), u3A.n23(), u3B.n13(), u3B.n23(),
        u3C.n13(), u3C.n23(), u3D.n13(), u3D.n23(),
        u3AB.n13(), u3AB.n23(), u3CD.n13(), u3CD.n23(),
        u3AC.n13(), u3AC.n23(), u3BD.n13(), u3BD.n23(),
        u3ABCD.n13(), u3ABCD.n23()};

      uint64_t key_1 = 
        (uint64_t(safe_key[0]) << 56) | (uint64_t(safe_key[1]) << 48) |
        (uint64_t(safe_key[2]) << 40) | (uint64_t(safe_key[3]) << 32) |
        (uint64_t(safe_key[4]) << 24) | (uint64_t(safe_key[5]) << 16) |
        (uint64_t(safe_key[6]) << 8) | uint64_t(safe_key[7]);

      uint64_t key_2 = 
        (uint64_t(safe_key[8]) << 56) | (uint64_t(safe_key[9]) << 48) |
        (uint64_t(safe_key[10]) << 40) | (uint64_t(safe_key[11]) << 32) |
        (uint64_t(safe_key[12]) << 24) | (uint64_t(safe_key[13]) << 16) |
        (uint64_t(safe_key[14]) << 8) | uint64_t(safe_key[15]);

      uint16_t key_3 = 
        (uint16_t(safe_key[16]) << 8) | uint16_t(safe_key[17]);

      key_type key = {key_1, key_2, key_3};
      /*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5],
                      safe_key[6], safe_key[7], safe_key[8],
                      safe_key[9], safe_key[10], safe_key[11],
                      safe_key[12], safe_key[13], safe_key[14],
                      safe_key[15], safe_key[16], safe_key[17]};
      */
      auto producer = [&u3A, &u3B, &u3C, &u3D, &u3AB, &u3CD, &u3AC, &u3BD, &u3ABCD, 
                       rhoABmax, rhoCDmax, rhoABCDmax,
                       rhoACmax, rhoBDmax, rhoACBDmax](key_type)
      { 
        return compute(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD, 
                       rhoABmax, rhoCDmax, rhoABCDmax,
                       rhoACmax, rhoBDmax, rhoACBDmax);
      };

      auto accessor = [](const auto& element) { *p_result_ = element.second; };

      cache_.visit(key, producer, accessor);
#else
      *p_result_ = compute(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD, 
                           rhoABmax, rhoCDmax, rhoABCDmax,
                           rhoACmax, rhoBDmax, rhoACBDmax);
#endif

      return *p_result_;
    }
    
    friend U3U9_db& u3u9db();

#ifdef U3LIB_ENABLE_U3U9_CACHE
    void clear() { cache_.clear(); }
#endif
  };

  inline U3U9_db& u3u9db()
  {
    static U3U9_db u3u9db;

    return u3u9db;
  }

} // namespace u3

#endif
