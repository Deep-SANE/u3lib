#ifndef UTILS_SO3_FANO_H
#define UTILS_SO3_FANO_H

#define UTILS_SO3_FANO_CACHE

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <utility>

#include <Int.h>
#include <so3Coeffs/SO3_CGC.h>          // namespace u2 (selection rule)
#include <so3Coeffs/wigner_backend.h>   // utils::wigner::wig9jj

#ifdef UTILS_SO3_FANO_CACHE
//#include <boost/unordered/unordered_flat_map.hpp>
#include <unordered_map>
#include <boost/container_hash/hash.hpp>
#endif

#ifdef UTILS_SO3_LOG
#include <fstream>
#endif

namespace utils
{
  // Fano (Racah-9j) coefficient: recoupling of 4 angular momenta.
  class SO3_Fano // recoupling of 4 angular momentum
  {
    // using key_type = std::array<uint8_t, 9>;
    using key_type = std::pair<uint64_t, uint8_t>;
    using value_type = double;

    enum { LL1, LL2, LL12, LL3, LL4, LL34, LL13, LL24, LL };

#ifdef UTILS_SO3_FANO_CACHE
    /*
    using map_type =
        boost::unordered_flat_map<key_type, value_type, boost::hash<key_type>>;
    */
    using map_type =
        std::unordered_map<key_type, value_type, boost::hash<key_type>>;

    inline static map_type* cache_ = nullptr;
#pragma omp threadprivate(cache_)
#endif

    static int to_int(auto i) { return static_cast<int>(i); }

    static value_type compute(const std::array<utils::UInt8, 9>& key)
    {
      int ll1 = to_int(key[LL1]);
      int ll2 = to_int(key[LL2]);
      int ll12 = to_int(key[LL12]);

      int ll3 = to_int(key[LL3]);
      int ll4 = to_int(key[LL4]);
      int ll34 = to_int(key[LL34]);

      int ll13 = to_int(key[LL13]);
      int ll24 = to_int(key[LL24]);
      int ll = to_int(key[LL]);

      double factor = std::sqrt((static_cast<double>(ll12) + 1.0) *
                                (static_cast<double>(ll34) + 1.0) *
                                (static_cast<double>(ll13) + 1.0) *
                                (static_cast<double>(ll24) + 1.0));

      double val9j = factor * wigner::wig9jj(ll1, ll2, ll12, ll3, ll4, ll34,
                                              ll13, ll24, ll);

      return val9j;
    }

    // This class is a singleton, and cannot be instantiated directly:
    SO3_Fano() = default;
    SO3_Fano(const SO3_Fano&) = delete;
    SO3_Fano& operator=(const SO3_Fano&) = delete;

  public:
    value_type get(std::array<utils::UInt8, 9> safe_key)
    {
#ifdef UTILS_SO3_LOG
      static std::ofstream f_log("so3_fano.log");

      for (auto i : safe_key)
        f_log << i << " ";
      f_log << "\n";
#endif

      // SU(2) selection rule:
      if ((u2::mult(safe_key[LL1], safe_key[LL2], safe_key[LL12]) == 0) ||
          (u2::mult(safe_key[LL3], safe_key[LL4], safe_key[LL34]) == 0) ||
          (u2::mult(safe_key[LL1], safe_key[LL3], safe_key[LL13]) == 0) ||
          (u2::mult(safe_key[LL2], safe_key[LL4], safe_key[LL24]) == 0) ||
          (u2::mult(safe_key[LL12], safe_key[LL34], safe_key[LL]) == 0) ||
          (u2::mult(safe_key[LL13], safe_key[LL24], safe_key[LL]) == 0))
        return 0.0;

      uint64_t key1 =
        (uint64_t(safe_key[0]) << 56) | (uint64_t(safe_key[1]) << 48) |
        (uint64_t(safe_key[2]) << 40) | (uint64_t(safe_key[3]) << 32) |
        (uint64_t(safe_key[4]) << 24) | (uint64_t(safe_key[5]) << 16) |
        (uint64_t(safe_key[6]) << 8) | uint64_t(safe_key[7]);

      key_type key = {key1, safe_key[8]};
/*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5],
                      safe_key[6], safe_key[7], safe_key[8]};
*/

#ifdef UTILS_SO3_FANO_CACHE
      auto iter = cache_->find(key);
      if (iter == cache_->end())
        iter = cache_->emplace(key, compute(safe_key)).first;
      return iter->second;
#else
      return compute(safe_key);
#endif
    }

#ifdef UTILS_SO3_FANO_CACHE
    static void initialize_cache()
    {
      assert(cache_ == nullptr);
      cache_ = new map_type;
    }

    static void finalize_cache()
    {
      assert(cache_ != nullptr);
      delete cache_;
      cache_ = nullptr;
    }
#endif

    friend SO3_Fano& fanos();
  };

  inline SO3_Fano& fanos()
  {
    static SO3_Fano fanos;

    return fanos;
  }

} // namespace utils

#endif
