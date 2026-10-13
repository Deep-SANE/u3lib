#ifndef UTILS_SO3_CGC_H
#define UTILS_SO3_CGC_H

#define UTILS_SO3_CGS_CACHE

#include <array>
#include <cmath>
#include <cstdint>

#include <Int.h>
#include <so3Coeffs/wigner_backend.h>

#ifdef UTILS_SO3_CGS_CACHE
//#include <boost/unordered/unordered_flat_map.hpp>
#include <unordered_map>
#endif

#ifdef UTILS_SO3_LOG
#include <fstream>
#endif

namespace u2
{
  inline int mult(int ll1, int ll2, int ll3)
  {
    if (std::abs(ll1 - ll2) <= ll3 && ll3 <= (ll1 + ll2))
      return 1;
    else
      return 0;
  }
} // namespace u2

namespace utils
{
  class SO3_CGC
  {
    // using key_type = std::array<int8_t, 6>;
    using key_type = uint64_t;
    using value_type = double;

    enum { LL1, LL2, LL3, MM1, MM2, MM3 };

#ifdef UTILS_SO3_CGS_CACHE
    // using map_type = boost::unordered_flat_map<key_type, value_type>;
    using map_type = std::unordered_map<key_type, value_type>;

    inline static map_type* cache_ = nullptr;
#pragma omp threadprivate(cache_)
#endif

    static int to_int(auto i) { return static_cast<int>(i); }

    static value_type compute(const std::array<utils::Int8, 6>& key)
    {
      int ll1 = to_int(key[LL1]);
      int ll2 = to_int(key[LL2]);
      int ll3 = to_int(key[LL3]);

      int mm1 = to_int(key[MM1]);
      int mm2 = to_int(key[MM2]);
      int mm3 = to_int(key[MM3]);

      double phase = ((-ll1 + ll2 - mm3) % 4 == 0) ? 1.0 : -1.0;
      double factor = std::sqrt(static_cast<double>(ll3) + 1.0);

      double val3j = phase * factor *
                     utils::wigner::wig3jj(ll1, ll2, ll3, mm1, mm2, -mm3);

      return val3j;
    }

    // This class is a singleton, and cannot be instantiated directly:
    SO3_CGC() = default;
    SO3_CGC(const SO3_CGC&) = delete;
    SO3_CGC& operator=(const SO3_CGC&) = delete;

  public:
    value_type get(std::array<utils::Int8, 6> safe_key)
    {
#ifdef UTILS_SO3_LOG
      static std::ofstream f_log("so3_cgc.log");

      for (auto i : safe_key)
        f_log << i << " ";
      f_log << "\n";
#endif

      // SU(2) selection rule:
      if ((u2::mult((int8_t)safe_key[LL1], (int8_t)safe_key[LL2],
                    (int8_t)safe_key[LL3]) == 0) ||
          (((int8_t)safe_key[MM1] + (int8_t)safe_key[MM2]) !=
           (int8_t)safe_key[MM3]))
        return 0.0;

      key_type key =
        (uint64_t(uint8_t(safe_key[0])) << 56) |
        (uint64_t(uint8_t(safe_key[1])) << 48) |
        (uint64_t(uint8_t(safe_key[2])) << 40) |
        (uint64_t(uint8_t(safe_key[3])) << 32) |
        (uint64_t(uint8_t(safe_key[4])) << 24) |
        (uint64_t(uint8_t(safe_key[5])) << 16) | 0;
/*
      key_type key = {safe_key[0], safe_key[1], safe_key[2],
                      safe_key[3], safe_key[4], safe_key[5]};
*/

#ifdef UTILS_SO3_CGS_CACHE
      auto iter = cache_->find(key);
      if (iter == cache_->end())
        iter = cache_->emplace(key, compute(safe_key)).first;
      return iter->second;
#else
      return compute(safe_key);
#endif
    }

#ifdef UTILS_SO3_CGS_CACHE
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

    friend SO3_CGC& cgcs();
  };

  inline SO3_CGC& cgcs()
  {
    static SO3_CGC cgcs;

    return cgcs;
  }

} // namespace utils

#endif
