#ifndef U3LIB_U3_U6Z6COEFFS_H
#define U3LIB_U3_U6Z6COEFFS_H

#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#ifdef TIME_U3U6Z6
#include <chrono>
#endif

#include <Eigen/Dense>

#include "U3basis.h"
#include "U3U2xU3U2.h"

#include <qcx_Cache.h>
#include <qcx_Int.h>
#include <so3Coeffs/SO3_CGC.h>
#include <so3Coeffs/SO3_Racah.h>

namespace u3
{
  // Input by users: u3A, u3B, u3ABC, u3C, u3AB, u3BC, all normalized
  inline auto
  U6(const u3::irrep& u3A, const std::vector<u3::gelfand>& gArrA,
     const u3::irrep& u3B, const std::vector<u3::gelfand>& gArrB,
     const u3::irrep& u3ABC, const std::vector<u3::gelfand>& gArrABC,
     const u3::irrep& u3C, const std::vector<u3::gelfand>& gArrC,
     const u3::irrep& u3AB, const std::vector<u3::gelfand>& gArrAB,
     int rhoABmax, int rhoAB_Cmax, const u3::irrep& u3BC,
     const std::vector<u3::gelfand>& gArrBC, int rhoBCmax, int rhoA_BCmax)
  {
#ifdef TIME_U3U6Z6
    static std::chrono::duration<double> u3u6_time(0);
    static double u3u6_last_print = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
#endif

    auto& racahs = qcx_utils::racahs();
    auto& wigU3U2db = u3::wigU3U2db();

    Eigen::MatrixXd A(rhoA_BCmax, rhoA_BCmax);
    Eigen::MatrixXd B(rhoA_BCmax, rhoABmax * rhoAB_Cmax * rhoBCmax);
    B.setZero();

    int r = 0;
    bool done = false;
    for (int n12A = u3A.n13(); n12A >= u3A.n23() && !done; n12A--)
    {
      for (int n22A = u3A.n23(); n22A >= u3A.n33() && !done; n22A--)
      {
        auto accA_BC = wigU3U2db.get(
            u3A.n13(), u3A.n23(), u3A.n33(), n12A, n22A,
            u3BC.n13(), u3BC.n23(), u3BC.n33(), u3BC.n13(), u3BC.n23(),
            u3ABC.n13(), u3ABC.n23(), u3ABC.n33(), u3ABC.n13(), u3ABC.n23());
        if (accA_BC.size() == 0) continue;
        
        A.row(r) = accA_BC;

        for (int n12B = u3B.n13(); n12B >= u3B.n23(); n12B--)
        {
          for (int n22B = u3B.n23(); n22B >= u3B.n33(); n22B--)
          {
            for (int n12AB = u3AB.n13(); n12AB >= u3AB.n23(); n12AB--)
            {
              for (int n22AB = u3AB.n23(); n22AB >= u3AB.n33(); n22AB--)
              {
                auto accAB = wigU3U2db.get(
                  u3A.n13(), u3A.n23(), u3A.n33(), n12A, n22A,
                  u3B.n13(), u3B.n23(), u3B.n33(), n12B, n22B,
                  u3AB.n13(), u3AB.n23(), u3AB.n33(), n12AB, n22AB);
                if (accAB.size() == 0) continue;

                for (int n12C = u3C.n13(); n12C >= u3C.n23(); n12C--)
                {
                  for (int n22C = u3C.n23(); n22C >= u3C.n33(); n22C--)
                  {
                    auto accAB_C = wigU3U2db.get(
                      u3AB.n13(), u3AB.n23(), u3AB.n33(), n12AB, n22AB,
                      u3C.n13(), u3C.n23(), u3C.n33(), n12C, n22C,
                      u3ABC.n13(), u3ABC.n23(), u3ABC.n33(), u3ABC.n13(), u3ABC.n23());
                    if (accAB_C.size() == 0) continue;

                    auto accBC = wigU3U2db.get(
                      u3B.n13(), u3B.n23(), u3B.n33(), n12B, n22B,
                      u3C.n13(), u3C.n23(), u3C.n33(), n12C, n22C,
                      u3BC.n13(), u3BC.n23(), u3BC.n33(), u3BC.n13(), u3BC.n23());
                    if (accBC.size() == 0) continue;

                    int llA = n12A - n22A, ll = u3ABC.n13() - u3ABC.n23(), 
                        llBC = u3BC.n13() - u3BC.n23(), llB = n12B - n22B, 
                        llAB = n12AB - n22AB, llC = n12C - n22C;
                    double racah = racahs.get({llA, llB, ll, llC, llAB, llBC});

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
        r++;

        if (r == rhoA_BCmax) done = true;
      }
    }

    // std::cout << A << std::endl;
    // std::cout << B << std::endl;

    Eigen::MatrixXd X = A.householderQr().solve(B);

    // std::vector<double> U(rhoABmax * rhoAB_Cmax * rhoBCmax * rhoA_BCmax);
    U3U6<coeffs_t> U(rhoABmax, rhoAB_Cmax, rhoBCmax, rhoA_BCmax);

    for (int rAB = 0; rAB < rhoABmax; rAB++)
      for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
        for (int rBC = 0; rBC < rhoBCmax; rBC++)
        {
          int col = (rAB * rhoAB_Cmax + rAB_C) * rhoBCmax + rBC;

          for (int rA_BC = 0; rA_BC < rhoA_BCmax; rA_BC++)
          {
            /*
            if (std::fabs(X(rA_BC, col)) < 1e-14)
              X(rA_BC, col) = 0.0;

            U[((rAB * rhoAB_Cmax + rAB_C) * rhoBCmax + rBC) * rhoA_BCmax +
              rA_BC] = X(rA_BC, col);
            */
            U[col * rhoA_BCmax + rA_BC] = 
              (std::fabs(X(rA_BC, col)) < 1e-14) ? 0.0 : X(rA_BC, col);
          }
        }

#ifdef TIME_U3U6Z6
    auto end = std::chrono::high_resolution_clock::now();
    u3u6_time += end - start;  // duration in seconds

    // Print occasionally 
    if (u3u6_time.count() - u3u6_last_print >= 10.0) 
    {
        std::cout << "Accumulated U3U6 time: " << u3u6_time.count() << " s\n";
        u3u6_last_print = u3u6_time.count();
    }
#endif

    return U;
  }

  // Database of U6 coefficients
  class U3U6_db
  {
    // using key_type = std::array<uint8_t, 12>;
    using key_type = std::pair<uint64_t, uint32_t>;
    using value_type = coeffs_t;

#ifdef U3LIB_ENABLE_U3U6_CACHE
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
    U3U6_db() = default;
    U3U6_db(const U3U6_db&) = delete;
    U3U6_db& operator=(const U3U6_db&) = delete;

    static value_type compute(
      const u3::irrep& u3A, const u3::irrep& u3B, 
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
                u3BC,gArrBC,rhoBCmax,rhoA_BCmax).data();
    }

  public:
    // pass by references to avoid copies
    U3U6<value_type> get_U3U6(
      const u3::irrep& u3A, const u3::irrep& u3B, 
      const u3::irrep& u3ABC, const u3::irrep& u3C,
      const u3::irrep& u3AB, int rhoABmax, int rhoAB_Cmax,
      const u3::irrep& u3BC, int rhoBCmax, int rhoA_BCmax)
    {
      auto coeffs = get(u3A, u3B, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                              u3BC, rhoBCmax, rhoA_BCmax);

      return {std::move(coeffs),
              static_cast<size_t>(rhoABmax),
              static_cast<size_t>(rhoAB_Cmax),
              static_cast<size_t>(rhoBCmax),
              static_cast<size_t>(rhoA_BCmax)};
    }

    // pass by values to normalize irreps
    value_type get(u3::irrep u3A, u3::irrep u3B, u3::irrep u3ABC, u3::irrep u3C, 
                    u3::irrep u3AB, int rhoABmax, int rhoAB_Cmax,
                    u3::irrep u3BC, int rhoBCmax, int rhoA_BCmax)
    {
      if (rhoABmax == 0 || rhoAB_Cmax == 0 || rhoBCmax == 0 || rhoA_BCmax == 0)
        return {};

      u3A.normalize();
      u3B.normalize();
      u3C.normalize();
      u3AB.normalize();
      u3BC.normalize();
      u3ABC.normalize();

      assert(rhoABmax == Rho(u3A,u3B,u3AB)); 
      assert(rhoAB_Cmax == Rho(u3AB,u3C,u3ABC));
      assert(rhoBCmax == Rho(u3B,u3C,u3BC)); 
      assert(rhoA_BCmax == Rho(u3A,u3BC,u3ABC));
      
#ifdef U3LIB_ENABLE_U3U6_CACHE
      std::array<qcx_utils::UInt8, 12> safe_key = {
        u3A.n13(), u3A.n23(), u3B.n13(), u3B.n23(),
        u3ABC.n13(), u3ABC.n23(), u3C.n13(), u3C.n23(),
        u3AB.n13(), u3AB.n23(), u3BC.n13(), u3BC.n23()};

      uint64_t key_1 = 
        (uint64_t(safe_key[0]) << 56) | (uint64_t(safe_key[1]) << 48) |
        (uint64_t(safe_key[2]) << 40) | (uint64_t(safe_key[3]) << 32) |
        (uint64_t(safe_key[4]) << 24) | (uint64_t(safe_key[5]) << 16) |
        (uint64_t(safe_key[6]) << 8) | uint64_t(safe_key[7]);

      uint32_t key_2 = 
        (uint32_t(safe_key[8]) << 24) | (uint32_t(safe_key[9]) << 16) |
        (uint32_t(safe_key[10]) << 8) | uint32_t(safe_key[11]);

      key_type key = {key_1, key_2};
      /*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5],
                      safe_key[6], safe_key[7], safe_key[8],
                      safe_key[9], safe_key[10], safe_key[11]};
      */

      auto producer = [&u3A, &u3B, &u3ABC, &u3C, 
                       &u3AB, &rhoABmax, &rhoAB_Cmax, 
                       &u3BC, &rhoBCmax, &rhoA_BCmax](key_type)
      { 
        return compute(u3A, u3B, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                             u3BC, rhoBCmax, rhoA_BCmax);
      };

      auto accessor = [](const auto& element) { *p_result_ = element.second; };

      cache_.visit(key, producer, accessor);
#else
      *p_result_ = compute(u3A, u3B, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                                 u3BC, rhoBCmax, rhoA_BCmax);
#endif

      return *p_result_;
    }

    // TO BE REMOVED IN THE FUTURE:
    value_type get(
      int n13a, int n23a, int n33a, int n13b, int n23b, int n33b, 
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

      if (rhoABmax == 0 || rhoAB_Cmax == 0 || rhoBCmax == 0 || rhoA_BCmax == 0)
        return {}; 
      
#ifdef U3LIB_ENABLE_U3U6_CACHE
      std::array<qcx_utils::UInt8, 12> safe_key = {
        u3A.n13(), u3A.n23(), u3B.n13(), u3B.n23(),
        u3ABC.n13(), u3ABC.n23(), u3C.n13(), u3C.n23(),
        u3AB.n13(), u3AB.n23(), u3BC.n13(), u3BC.n23()};

      uint64_t key_1 = 
        (uint64_t(safe_key[0]) << 56) | (uint64_t(safe_key[1]) << 48) |
        (uint64_t(safe_key[2]) << 40) | (uint64_t(safe_key[3]) << 32) |
        (uint64_t(safe_key[4]) << 24) | (uint64_t(safe_key[5]) << 16) |
        (uint64_t(safe_key[6]) << 8) | uint64_t(safe_key[7]);

      uint32_t key_2 = 
        (uint32_t(safe_key[8]) << 24) | (uint32_t(safe_key[9]) << 16) |
        (uint32_t(safe_key[10]) << 8) | uint32_t(safe_key[11]);

      key_type key = {key_1, key_2};
      /*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5],
                      safe_key[6], safe_key[7], safe_key[8],
                      safe_key[9], safe_key[10], safe_key[11]};
      */

      auto producer = [&u3A, &u3B, &u3ABC, &u3C, 
                       &u3AB, &rhoABmax, &rhoAB_Cmax, 
                       &u3BC, &rhoBCmax, &rhoA_BCmax](key_type)
      { 
        return compute(u3A, u3B, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                             u3BC, rhoBCmax, rhoA_BCmax);
      };

      auto accessor = [](const auto& element) { *p_result_ = element.second; };

      cache_.visit(key, producer, accessor);
#else
      *p_result_ = compute(u3A, u3B, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                                 u3BC, rhoBCmax, rhoA_BCmax);
#endif

      return *p_result_;
    }
  
    friend U3U6_db& u3u6db();

#ifdef U3LIB_ENABLE_U3U6_CACHE
    void clear() { cache_.clear(); }
#endif
  };

  inline U3U6_db& u3u6db()
  {
    static U3U6_db u3u6db;

    return u3u6db; 
  }

  // Input by users: u3B, u3A, u3ABC, u3C, u3AB, u3AC
  inline auto
  Z6(const u3::irrep& u3B, const std::vector<u3::gelfand>& gArrB,
     const u3::irrep& u3A, const std::vector<u3::gelfand>& gArrA,
     const u3::irrep& u3ABC, const std::vector<u3::gelfand>& gArrABC,
     const u3::irrep& u3C, const std::vector<u3::gelfand>& gArrC,
     const u3::irrep& u3AB, const std::vector<u3::gelfand>& gArrAB,
     int rhoABmax, int rhoAB_Cmax, const u3::irrep& u3AC,
     const std::vector<u3::gelfand>& gArrAC, int rhoACmax, int rhoAC_Bmax)
  {
#ifdef TIME_U3U6Z6
    static std::chrono::duration<double> u3z6_time(0);
    static double u3z6_last_print = 0.0;
    auto start = std::chrono::high_resolution_clock::now();
#endif

    auto& racahs = qcx_utils::racahs();
    auto& wigU3U2db = u3::wigU3U2db();

    Eigen::MatrixXd A(rhoAC_Bmax, rhoAC_Bmax);
    Eigen::MatrixXd B(rhoAC_Bmax, rhoABmax * rhoAB_Cmax * rhoACmax);
    B.setZero();

    int r = 0;
    bool done = false;
    for (int n12B = u3B.n13(); n12B >= u3B.n23() && !done; n12B--)
    {
      for (int n22B = u3B.n23(); n22B >= u3B.n33() && !done; n22B--)
      {
        auto accAC_B = wigU3U2db.get(
          u3AC.n13(), u3AC.n23(), u3AC.n33(), u3AC.n13(), u3AC.n23(), 
          u3B.n13(), u3B.n23(), u3B.n33(), n12B, n22B,
          u3ABC.n13(), u3ABC.n23(), u3ABC.n33(), u3ABC.n13(), u3ABC.n23());
        if (accAC_B.size() == 0) continue;

        A.row(r) = accAC_B;

        for (int n12A = u3A.n13(); n12A >= u3A.n23(); n12A--)
        {
          for (int n22A = u3A.n23(); n22A >= u3A.n33(); n22A--)
          {
            for (int n12AB = u3AB.n13(); n12AB >= u3AB.n23(); n12AB--)
            {
              for (int n22AB = u3AB.n23(); n22AB >= u3AB.n33(); n22AB--)
              {
                auto accAB = wigU3U2db.get(
                  u3A.n13(), u3A.n23(), u3A.n33(), n12A, n22A, 
                  u3B.n13(), u3B.n23(), u3B.n33(), n12B, n22B, 
                  u3AB.n13(), u3AB.n23(), u3AB.n33(), n12AB, n22AB);
                if (accAB.size() == 0) continue;
                
                for (int n12C = u3C.n13(); n12C >= u3C.n23(); n12C--)
                {
                  for (int n22C = u3C.n23(); n22C >= u3C.n33(); n22C--)
                  {
                    auto accAB_C = wigU3U2db.get(
                      u3AB.n13(), u3AB.n23(), u3AB.n33(), n12AB, n22AB, 
                      u3C.n13(), u3C.n23(), u3C.n33(), n12C, n22C, 
                      u3ABC.n13(), u3ABC.n23(), u3ABC.n33(), u3ABC.n13(), u3ABC.n23());
                    if (accAB_C.size() == 0) continue;

                    auto accAC = wigU3U2db.get(
                      u3A.n13(), u3A.n23(), u3A.n33(), n12A, n22A, 
                      u3C.n13(), u3C.n23(), u3C.n33(), n12C, n22C,
                      u3AC.n13(), u3AC.n23(), u3AC.n33(), u3AC.n13(), u3AC.n23());
                    if (accAC.size() == 0) continue;

                    int llA = n12A - n22A, llB = n12B - n22B,
                        ll = u3ABC.n13() - u3ABC.n23(), llC = n12C - n22C,
                        llAB = n12AB - n22AB, llAC = u3AC.n13() - u3AC.n23();

                    double phase = ((llA + ll - llAB - llAC) % 4 == 0) ? 1 : -1;
                    double racah = racahs.get({llB, llA, ll, llC, llAB, llAC});

                    for (int rAB = 0; rAB < rhoABmax; rAB++)
                      for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
                        for (int rAC = 0; rAC < rhoACmax; rAC++)
                          B(r, (rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC) 
                            += phase * racah * accAB(rAB) * accAB_C(rAB_C) * accAC(rAC);
                  }
                }
              }
            }
          }
        }
        r++;

        if (r == rhoAC_Bmax) done = true;
      }
    }

    // std::cout << A << std::endl;
    // std::cout << B << std::endl;

    Eigen::MatrixXd X = A.householderQr().solve(B);

    // std::vector<double> Z(rhoABmax * rhoAB_Cmax * rhoACmax * rhoAC_Bmax);
    U3Z6<coeffs_t> Z(rhoABmax, rhoAB_Cmax, rhoACmax, rhoAC_Bmax);

    for (int rAB = 0; rAB < rhoABmax; rAB++)
      for (int rAB_C = 0; rAB_C < rhoAB_Cmax; rAB_C++)
        for (int rAC = 0; rAC < rhoACmax; rAC++)
        {
          int col = (rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC;
          for (int rAC_B = 0; rAC_B < rhoAC_Bmax; rAC_B++)
          {
            /*
            if (std::fabs(X(rAC_B, col)) < 1e-14)
              X(rAC_B, col) = 0.0;

            Z[((rAB * rhoAB_Cmax + rAB_C) * rhoACmax + rAC) * rhoAC_Bmax +
              rAC_B] = X(rAC_B, col);
            */
            Z[col * rhoAC_Bmax + rAC_B] =
              (std::fabs(X(rAC_B, col)) < 1e-14) ? 0.0 : X(rAC_B, col);
          }
        }

#ifdef TIME_U3U6Z6
    auto end = std::chrono::high_resolution_clock::now();
    u3z6_time += end - start;  // duration in seconds

    // Print occasionally 
    if (u3z6_time.count() - u3z6_last_print >= 10.0) 
    {
        std::cout << "Accumulated U3Z6 time: " << u3z6_time.count() << " s\n";
        u3z6_last_print = u3z6_time.count();
    }
#endif

    return Z;
  }

  // Database of Z6 coefficients
  class U3Z6_db
  {
    // using key_type = std::array<uint8_t, 12>;
    using key_type = std::pair<uint64_t, uint32_t>;
    using value_type = coeffs_t;

#ifdef U3LIB_ENABLE_U3Z6_CACHE
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
    U3Z6_db() = default;
    U3Z6_db(const U3Z6_db&) = delete;
    U3Z6_db& operator=(const U3Z6_db&) = delete;

    static value_type compute(
      const u3::irrep& u3B, const u3::irrep& u3A, 
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

      return Z6(u3B,gArrB,u3A,gArrA,u3ABC,gArrABC,u3C,gArrC,
                u3AB,gArrAB,rhoABmax,rhoAB_Cmax,
                u3AC,gArrAC,rhoACmax,rhoAC_Bmax).data();
    }

  public:
    // pass by references to avoid copies
    U3Z6<value_type> get_U3Z6(
      const u3::irrep& u3B, const u3::irrep& u3A, 
      const u3::irrep& u3ABC, const u3::irrep& u3C,
      const u3::irrep& u3AB, int rhoABmax, int rhoAB_Cmax,
      const u3::irrep& u3AC, int rhoACmax, int rhoAC_Bmax)
    {
      auto coeffs = get(u3B, u3A, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                              u3AC, rhoACmax, rhoAC_Bmax);

      return {std::move(coeffs),
              static_cast<size_t>(rhoABmax),
              static_cast<size_t>(rhoAB_Cmax),
              static_cast<size_t>(rhoACmax),
              static_cast<size_t>(rhoAC_Bmax)};
    }

    // pass by values to normalize irreps
    value_type get(u3::irrep u3B, u3::irrep u3A, u3::irrep u3ABC, u3::irrep u3C,
                   u3::irrep u3AB, int rhoABmax, int rhoAB_Cmax,
                   u3::irrep u3AC, int rhoACmax, int rhoAC_Bmax)
    {
      if (rhoABmax == 0 || rhoAB_Cmax == 0 || rhoACmax == 0 || rhoAC_Bmax == 0)
        return {};

      u3A.normalize();
      u3B.normalize();
      u3C.normalize();
      u3AB.normalize();
      u3AC.normalize();
      u3ABC.normalize();

      assert(rhoABmax == Rho(u3A,u3B,u3AB)); 
      assert(rhoAB_Cmax == Rho(u3AB,u3C,u3ABC));
      assert(rhoACmax == Rho(u3A,u3C,u3AC)); 
      assert(rhoAC_Bmax == Rho(u3AC,u3B,u3ABC));

#ifdef U3LIB_ENABLE_U3Z6_CACHE
      std::array<qcx_utils::UInt8, 12> safe_key = {
        u3B.n13(), u3B.n23(), u3A.n13(), u3A.n23(),
        u3ABC.n13(), u3ABC.n23(), u3C.n13(), u3C.n23(),
        u3AB.n13(), u3AB.n23(), u3AC.n13(), u3AC.n23()};

      uint64_t key_1 = 
        (uint64_t(safe_key[0]) << 56) | (uint64_t(safe_key[1]) << 48) |
        (uint64_t(safe_key[2]) << 40) | (uint64_t(safe_key[3]) << 32) |
        (uint64_t(safe_key[4]) << 24) | (uint64_t(safe_key[5]) << 16) |
        (uint64_t(safe_key[6]) << 8) | uint64_t(safe_key[7]);

      uint32_t key_2 = 
        (uint32_t(safe_key[8]) << 24) | (uint32_t(safe_key[9]) << 16) |
        (uint32_t(safe_key[10]) << 8) | uint32_t(safe_key[11]);

      key_type key = {key_1, key_2};
      /*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5],
                      safe_key[6], safe_key[7], safe_key[8],
                      safe_key[9], safe_key[10], safe_key[11]};
      */
      auto producer = [&u3B, &u3A, &u3ABC, &u3C, 
                       &u3AB, &rhoABmax, &rhoAB_Cmax,
                       &u3AC, &rhoACmax, &rhoAC_Bmax](key_type)
      { 
        return compute(u3B, u3A, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                             u3AC, rhoACmax, rhoAC_Bmax);
      };

      auto accessor = [](const auto& element) { *p_result_ = element.second; };

      cache_.visit(key, producer, accessor);
#else
      *p_result_ = compute(u3B, u3A, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                                 u3AC, rhoACmax, rhoAC_Bmax);
#endif

      return *p_result_;
    }

    // TO BE REMOVED IN THE FUTURE:
    value_type get(
      int n13b, int n23b, int n33b, int n13a, int n23a, int n33a,
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

      if (rhoABmax == 0 || rhoAB_Cmax == 0 || rhoACmax == 0 || rhoAC_Bmax == 0)
        return {}; 
      
#ifdef U3LIB_ENABLE_U3Z6_CACHE
      std::array<qcx_utils::UInt8, 12> safe_key = {
        u3B.n13(), u3B.n23(), u3A.n13(), u3A.n23(),
        u3ABC.n13(), u3ABC.n23(), u3C.n13(), u3C.n23(),
        u3AB.n13(), u3AB.n23(), u3AC.n13(), u3AC.n23()};

      uint64_t key_1 = 
        (uint64_t(safe_key[0]) << 56) | (uint64_t(safe_key[1]) << 48) |
        (uint64_t(safe_key[2]) << 40) | (uint64_t(safe_key[3]) << 32) |
        (uint64_t(safe_key[4]) << 24) | (uint64_t(safe_key[5]) << 16) |
        (uint64_t(safe_key[6]) << 8) | uint64_t(safe_key[7]);

      uint32_t key_2 = 
        (uint32_t(safe_key[8]) << 24) | (uint32_t(safe_key[9]) << 16) |
        (uint32_t(safe_key[10]) << 8) | uint32_t(safe_key[11]);

      key_type key = {key_1, key_2};
      /*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5],
                      safe_key[6], safe_key[7], safe_key[8],
                      safe_key[9], safe_key[10], safe_key[11]};
      */
      auto producer = [&u3B, &u3A, &u3ABC, &u3C, 
                       &u3AB, &rhoABmax, &rhoAB_Cmax,
                       &u3AC, &rhoACmax, &rhoAC_Bmax](key_type)
      { 
        return compute(u3B, u3A, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                             u3AC, rhoACmax, rhoAC_Bmax);
      };

      auto accessor = [](const auto& element) { *p_result_ = element.second; };

      cache_.visit(key, producer, accessor);
#else
      *p_result_ = compute(u3B, u3A, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                                 u3AC, rhoACmax, rhoAC_Bmax);
#endif

      return *p_result_;
    }

    friend U3Z6_db& u3z6db();

#ifdef U3LIB_ENABLE_U3Z6_CACHE
    void clear() { cache_.clear(); }
#endif
  };

  inline U3Z6_db& u3z6db()
  {
    static U3Z6_db u3z6db;

    return u3z6db;
  }

} // namespace u3
#endif
