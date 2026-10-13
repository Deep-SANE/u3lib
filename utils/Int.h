#ifndef UTILS_INT_H
#define UTILS_INT_H

#include <cassert>
#include <cstdint>
#include <concepts>
#include <iostream>
#include <limits>
#include <type_traits>

namespace utils
{
  namespace detail
  {
    template <std::integral I>
    constexpr int64_t as_int64_t(I i)
      requires(sizeof(I) <= sizeof(int64_t))
    {
      if constexpr (sizeof(I) < sizeof(int64_t))
        return static_cast<int64_t>(i);
      else
      {
        if constexpr (std::is_signed_v<I>)
          return i;
        else
        {
          int64_t ub =
              static_cast<uint64_t>(std::numeric_limits<int64_t>::max());

          assert(i <= ub);

          return static_cast<int64_t>(i);
        }
      }
    }

    template <std::integral I>
      requires(sizeof(I) < sizeof(int64_t))
    inline constexpr int64_t lbi64 = as_int64_t(std::numeric_limits<I>::min());

    template <std::integral I>
      requires(sizeof(I) < sizeof(int64_t))
    inline constexpr int64_t ubi64 = as_int64_t(std::numeric_limits<I>::max());

    int8_t to_int8_t(std::integral auto i)
    // requires(sizeof(i) < sizeof(int64_t))
    {
      assert((as_int64_t(i) >= lbi64<int8_t>) &&
             (as_int64_t(i) <= ubi64<int8_t>));

      return static_cast<int8_t>(i);
    }

    uint8_t to_uint8_t(std::integral auto i)
    //  requires(sizeof(i) < sizeof(int64_t))
    {
      assert((as_int64_t(i) >= lbi64<uint8_t>) &&
             (as_int64_t(i) <= ubi64<uint8_t>));

      return static_cast<uint8_t>(i);
    }
  } // namespace int

  struct Int8
  {
    int8_t value;

    Int8(std::integral auto i) : value(detail::to_int8_t(i)) { }

    operator int8_t() const { return value; }
  };

  inline std::ostream& operator<<(std::ostream& s, const Int8& x)
  {
    return s << static_cast<int>(x.value);
  }

  struct UInt8
  {
    uint8_t value;

    UInt8(std::integral auto i) : value(detail::to_uint8_t(i)) { }

    operator uint8_t() const { return value; }
  };

  inline std::ostream& operator<<(std::ostream& s, const UInt8& x)
  {
    return s << static_cast<int>(x.value);
  }

} // namespace utils

#endif
