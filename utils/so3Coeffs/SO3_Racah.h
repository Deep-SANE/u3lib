#ifndef UTILS_SO3_RACAH_H
#define UTILS_SO3_RACAH_H

#define UTILS_SO3_RACAH_CACHE

#include <array>
#include <cmath>
#include <cstdint>

#include <Int.h>
#include <so3Coeffs/SO3_CGC.h>
#include <so3Coeffs/wigner_backend.h>

#ifdef UTILS_SO3_RACAH_CACHE
//#include <boost/unordered/unordered_flat_map.hpp>
#include <unordered_map>
#endif

#ifdef UTILS_SO3_LOG
#include <fstream>
#endif

namespace utils
{
  class SO3_Racah // recoupling of 3 angular momentum
  {
    // using key_type = std::array<uint8_t, 6>;
    using key_type = uint64_t;
    using value_type = double;

    enum { LL1, LL2, LL, LL3, LL12, LL23 };

#ifdef UTILS_SO3_RACAH_CACHE
    // using map_type = boost::unordered_flat_map<key_type, value_type>;
    using map_type = std::unordered_map<key_type, value_type>;

    inline static map_type* cache_ = nullptr;
#pragma omp threadprivate(cache_)
#endif

    static int to_int(auto i) { return static_cast<int>(i); }

    static value_type compute(const std::array<utils::UInt8, 6>& key)
    {
      int ll1 = to_int(key[LL1]);
      int ll2 = to_int(key[LL2]);
      int ll = to_int(key[LL]);

      int ll3 = to_int(key[LL3]);
      int ll12 = to_int(key[LL12]);
      int ll23 = to_int(key[LL23]);

      double phase = ((ll1 + ll2 + ll3 + ll) % 4 == 0) ? 1 : -1;

      double factor = std::sqrt((static_cast<double>(ll12) + 1.0) *
                                (static_cast<double>(ll23) + 1.0));

      double val6j = phase * factor *
                     utils::wigner::wig6jj(ll1, ll2, ll12, ll3, ll, ll23);

      return val6j;
    }

    // This class is a singleton, and cannot be instantiated directly:
    SO3_Racah() = default;
    SO3_Racah(const SO3_Racah&) = delete;
    SO3_Racah& operator=(const SO3_Racah&) = delete;

  public:
    value_type get(std::array<utils::UInt8, 6> safe_key)
    {
#ifdef UTILS_SO3_LOG
      static std::ofstream f_log("so3_racah.log");

      for (auto i : safe_key)
        f_log << i << " ";
      f_log << "\n";
#endif

      // SU(2) selection rule:
      if ((u2::mult(safe_key[LL1], safe_key[LL2], safe_key[LL12]) == 0) || 
          (u2::mult(safe_key[LL12], safe_key[LL3], safe_key[LL]) == 0) ||
          (u2::mult(safe_key[LL2], safe_key[LL3], safe_key[LL23]) == 0) || 
          (u2::mult(safe_key[LL1], safe_key[LL23], safe_key[LL]) == 0))
        return 0.0;

      key_type key = 
        (uint64_t(safe_key[0]) << 56) | (uint64_t(safe_key[1]) << 48) |
        (uint64_t(safe_key[2]) << 40) | (uint64_t(safe_key[3]) << 32) |
        (uint64_t(safe_key[4]) << 24) | (uint64_t(safe_key[5]) << 16) | 0;
/*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5]};
*/

#ifdef UTILS_SO3_RACAH_CACHE
      auto iter = cache_->find(key);
      if (iter == cache_->end())
        iter = cache_->emplace(key, compute(safe_key)).first;
      return iter->second;
#else
      return compute(safe_key);
#endif
    }

#ifdef UTILS_SO3_RACAH_CACHE
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

    friend SO3_Racah& racahs();
  };

  inline SO3_Racah& racahs()
  {
    static SO3_Racah racahs;

    return racahs;
  }

} // namespace utils

#endif
