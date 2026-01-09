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
     int rhoABmax, int rhoCDmax, int rhoAB_CDmax, int rhoACmax, int rhoBDmax,
     int rhoAC_BDmax, qcx_utils::SO3_Racah& racahs,
     u3::WigU3U2_db& wigU3U2db)
  {
    auto& cgcs = qcx_utils::cgcs();

    std::vector<double> U(rhoACmax * rhoBDmax * rhoAB_CDmax * rhoABmax *
                              rhoCDmax * rhoAC_BDmax, 0.0);

    // couple u3AB and u3C
    auto [u3AB_C, rhomaxsAB_C] = u3::decomposition(u3AB, u3C);

    /* an array to store the couple I-II-III irreps
    std::vector<u3::irrep> u3ABC;
    std::vector<int> rhoAB_Cmax, rhoAC_Bmax;

    for (int i = 0; i < u3AB_C.size(); i++)
    {
      int rhomaxAC_B = u3::Rho(u3AC, u3B, u3AB_C[i]);
      if (rhomaxAC_B > 0)
      {
        u3ABC.push_back(u3AB_C[i]);
        rhoAB_Cmax.push_back(rhomaxsAB_C[i]);
        rhoAC_Bmax.push_back(rhomaxAC_B);
      }
    } */
    
    //for (int i = 0; i < u3ABC.size(); i++)
    for (int i = 0; i < u3AB_C.size(); i++)
    {
      int rhomaxAC_B = u3::Rho(u3AC, u3B, u3AB_C[i]);
      if (rhomaxAC_B == 0)
        continue;

      int rhoABC_Dmax = u3::Rho(u3AB_C[i], u3D, u3ABCD);
      if (rhoABC_Dmax == 0)
        continue;

      std::vector<u3::gelfand> gArrABC = u3::gelfandBasis(u3AB_C[i]);

      std::vector<double> U6_1 =
          u3::U6(u3AC, gArrAC, u3B, gArrB, u3ABCD, gArrABCD, u3D, gArrD,
                 u3AB_C[i], gArrABC, rhomaxAC_B, rhoABC_Dmax, u3BD, gArrBD,
                 rhoBDmax, rhoAC_BDmax, racahs, wigU3U2db);

      std::vector<double> U6_2 =
          u3::U6(u3AB, gArrAB, u3C, gArrC, u3ABCD, gArrABCD, u3D, gArrD,
                 u3AB_C[i], gArrABC, rhomaxsAB_C[i], rhoABC_Dmax, u3CD, gArrCD,
                 rhoCDmax, rhoAB_CDmax, racahs, wigU3U2db);

      std::vector<double> Z6 =
          u3::Z6(u3B, gArrB, u3A, gArrA, u3AB_C[i], gArrABC, u3C, gArrC, u3AB,
                 gArrAB, rhoABmax, rhomaxsAB_C[i], u3AC, gArrAC, rhoACmax,
                 rhomaxAC_B, racahs, wigU3U2db);

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
                  for (int rAC_B = 0; rAC_B < rhomaxAC_B; rAC_B++)
                  {
                    for (int rAB_C = 0; rAB_C < rhomaxsAB_C[i]; rAB_C++)
                    {
                      for (int rABC_D = 0; rABC_D < rhoABC_Dmax; rABC_D++)
                      {
                        U[((((rAC * rhoBDmax + rBD) * rhoAB_CDmax + rAB_CD) *
                        rhoABmax + rAB) * rhoCDmax + rCD) * rhoAC_BDmax + rAC_BD] 
                        += U6_1[((rAC_B * rhoABC_Dmax + rABC_D) * rhoBDmax +
                                  rBD) * rhoAC_BDmax + rAC_BD] 
                         * Z6[((rAB * rhomaxsAB_C[i] + rAB_C) * rhoACmax +
                                rAC) * rhomaxAC_B + rAC_B] 
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

  template <typename T>
  class U3U9
  {
    T coeffs_;
    qcx_utils::NDIndex6 ndi_;

  public:
    U3U9()
    {
      ndi_.add_dimensions(0, 0, 0, 0, 0, 0);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    U3U9(T&& coeffs, size_t rhoACmax, size_t rhoBDmax, size_t rhoABCDmax,
         size_t rhoABmax, size_t rhoCDmax, size_t rhoACBDmax)
        : coeffs_(std::move(coeffs))
    {
      ndi_.add_dimensions(rhoACmax, rhoBDmax, rhoABCDmax, rhoABmax, rhoCDmax,
                          rhoACBDmax);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    size_t size() const { return coeffs_.size(); }

    bool empty() const { return size() == 0; }

    size_t rhoACmax() const { return ndi_.dimension(0); }
    size_t rhoBDmax() const { return ndi_.dimension(1); }
    size_t rhoABCDmax() const { return ndi_.dimension(2); }
    size_t rhoABmax() const { return ndi_.dimension(3); }
    size_t rhoCDmax() const { return ndi_.dimension(4); }
    size_t rhoACBDmax() const { return ndi_.dimension(5); }

    double operator[](size_t i) const
    {
      assert(i < coeffs_.size());
      return coeffs_[i];
    }

    double get(size_t rhoAC, size_t rhoBD, size_t rhoABCD, size_t rhoAB,
               size_t rhoCD, size_t rhoACBD)
    {
      std::array multi_index = {rhoAC, rhoBD, rhoABCD, rhoAB, rhoCD, rhoACBD};
      size_t i = ndi_.to_linear(multi_index);
      return operator[](i);
    }
  };

#define U3LIB_ENABLE_U3U9_CACHE
  class U3U9_db
  {
#ifdef U3LIB_ENABLE_U3U9_CACHE
    using key_type = std::array<int, 18>;
    std::map<key_type, std::vector<double>> map_;

    using return_type = std::span<const double>;
#else
    using return_type = std::vector<double>;
#endif

    qcx_utils::SO3_CGC& cgcs_ = qcx_utils::cgcs();
    qcx_utils::SO3_Racah& racahs_;
    u3::WigU3U2_db& wigU3U2db_;

    auto calculate(const u3::irrep& u3A, const u3::irrep& u3B,
                   const u3::irrep& u3C, const u3::irrep& u3D,
                   const u3::irrep& u3AB, const u3::irrep& u3CD,
                   const u3::irrep& u3AC, const u3::irrep& u3BD,
                   const u3::irrep& u3ABCD, int rhoABmax, int rhoCDmax,
                   int rhoAB_CDmax, int rhoACmax, int rhoBDmax, int rhoAC_BDmax)
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
                rhoAC_BDmax, racahs_, wigU3U2db_);
    }

  public:
    U3U9_db(qcx_utils::SO3_Racah& racahs,
            u3::WigU3U2_db& wigU3U2db)
        : racahs_(racahs), wigU3U2db_(wigU3U2db)
    { }

    U3U9<return_type> get_U3U9(u3::irrep u3A, u3::irrep u3B, u3::irrep u3C,
                               u3::irrep u3D, u3::irrep u3AB, u3::irrep u3CD,
                               u3::irrep u3AC, u3::irrep u3BD, u3::irrep u3ABCD,
                               int rhoACmax, int rhoBDmax, int rhoABCDmax,
                               int rhoABmax, int rhoCDmax, int rhoACBDmax)
    {
      if ((rhoACmax == 0) || (rhoBDmax == 0) || (rhoABCDmax == 0) ||
          (rhoABmax == 0) || (rhoCDmax == 0) || (rhoACBDmax == 0))
        return {};

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

    /*
    auto get(u3::irrep u3A, u3::irrep u3B, u3::irrep u3C, u3::irrep u3D,
             u3::irrep u3AB, u3::irrep u3CD, u3::irrep u3AC, u3::irrep u3BD,
             u3::irrep u3ABCD)
    {
      auto rhoABmax = Rho(u3A, u3B, u3AB);
      auto rhoCDmax = Rho(u3C, u3D, u3CD);
      auto rhoABCDmax = Rho(u3AB, u3CD, u3ABCD);
      auto rhoACmax = Rho(u3A, u3C, u3AC);
      auto rhoBDmax = Rho(u3B, u3D, u3BD);
      auto rhoACBDmax = Rho(u3AC, u3BD, u3ABCD);

      return get(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD, rhoABmax,
                 rhoCDmax, rhoABCDmax, rhoACmax, rhoBDmax, rhoACBDmax);
    }
    */

    return_type get(u3::irrep u3A, u3::irrep u3B, u3::irrep u3C, u3::irrep u3D,
                    u3::irrep u3AB, u3::irrep u3CD, u3::irrep u3AC,
                    u3::irrep u3BD, u3::irrep u3ABCD, int rhoABmax,
                    int rhoCDmax, int rhoABCDmax, int rhoACmax, int rhoBDmax,
                    int rhoACBDmax)
    {
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

      if ((rhoABmax == 0) || (rhoCDmax == 0) || (rhoABCDmax == 0) ||
          (rhoACmax == 0) || (rhoBDmax == 0) || (rhoACBDmax == 0))
        return {};

#ifdef U3LIB_ENABLE_U3U9_CACHE
      key_type key = {u3A.n13(),    u3A.n23(),   u3B.n13(),  u3B.n23(),
                      u3C.n13(),    u3C.n23(),   u3D.n13(),  u3D.n23(),
                      u3AB.n13(),   u3AB.n23(),  u3CD.n13(), u3CD.n23(),
                      u3AC.n13(),   u3AC.n23(),  u3BD.n13(), u3BD.n23(),
                      u3ABCD.n13(), u3ABCD.n23()};

      auto iter = map_.find(key);
      if (iter == map_.end())
      {
        auto p = map_.try_emplace(key, calculate(u3A, u3B, u3C, u3D, u3AB, u3CD,
                                                 u3AC, u3BD, u3ABCD, rhoABmax,
                                                 rhoCDmax, rhoABCDmax, rhoACmax,
                                                 rhoBDmax, rhoACBDmax));

        assert(p.second);
        iter = p.first;
      }

      return iter->second;
#else
      return calculate(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD,
                       rhoABmax, rhoCDmax, rhoABCDmax, rhoACmax, rhoBDmax,
                       rhoACBDmax);
#endif
    }

    // TO BE REMOVED IN THE FUTURE:
    auto get(int n13a, int n23a, int n33a, int n13b, int n23b, int n33b, 
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
          rhoAB_CDmax = Rho(u3AB,u3CD,u3ABCD), 
          rhoAC_BDmax = Rho(u3AC,u3BD,u3ABCD);

      bool flag = false;
      std::vector<double> c;

      if (rhoABmax > 0 && rhoCDmax > 0 && rhoACmax > 0 && rhoBDmax > 0
          && rhoAB_CDmax > 0 && rhoAC_BDmax > 0)
      {
#ifdef U3LIB_ENABLE_U3U9_CACHE
        key_type key = {u3A.n13(), u3A.n23(),
                        u3B.n13(), u3B.n23(),
                        u3C.n13(), u3C.n23(),
                        u3D.n13(), u3D.n23(),
                        u3AB.n13(), u3AB.n23(),
                        u3CD.n13(), u3CD.n23(),
                        u3AC.n13(), u3AC.n23(),
                        u3BD.n13(), u3BD.n23(),
                        u3ABCD.n13(), u3ABCD.n23()};
        auto iter = map_.find(key);
        if (iter == map_.end())
        {
          auto p = map_.try_emplace(key, 
            calculate(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD, 
              rhoABmax,rhoCDmax,rhoAB_CDmax,rhoACmax,rhoBDmax,rhoAC_BDmax));
          assert(p.second);
          iter = p.first;
        }

        c = iter->second;
#else
        c = calculate(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD, 
              rhoABmax,rhoCDmax,rhoAB_CDmax,rhoACmax,rhoBDmax,rhoAC_BDmax);
#endif
        flag = true;
      }

      return std::make_pair(flag, c);
    }
  };

} // namespace u3

#endif
