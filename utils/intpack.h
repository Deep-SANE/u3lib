#ifndef UTILS_INTPACK_H
#define UTILS_INTPACK_H

#include <cassert>
#include <concepts>
#include <cstdint>
#include <cstring>
#include <limits>

#include <int128.h>

namespace utils
{
  namespace intpack
  {
    // Helper trait; the resulting unsigned type is allowed to have 128 bits
    template <typename T>
    concept unsigned_integral =
        // std::unsigned_integral<T> || std::same_as<T, utils::uint128_t>;
        std::unsigned_integral<T>;

    // Helper function to mask 'B' lowest bits for unsigned types:
    template <std::unsigned_integral U, unsigned B>
    constexpr U mask()
    {
      static_assert(B <= 8 * sizeof(U));

      if constexpr (B == 8 * sizeof(U))
        return std::numeric_limits<U>::max();
      else
        return (static_cast<U>(1) << B) - 1;
    }

    // Checks if the value of a signed integer fits into the lowest 'B' bits:
    template <unsigned B, std::signed_integral S>
    bool check(S s)
    {
      static_assert(B > 0);
      static_assert(B <= 8 * sizeof(S));

      static constexpr S lb = (B == sizeof(S) * 8)
                                  ? std::numeric_limits<S>::min()
                                  : -(static_cast<S>(1) << (B - 1));

      static constexpr S ub = (B == sizeof(S) * 8)
                                  ? std::numeric_limits<S>::max()
                                  : (static_cast<S>(1) << (B - 1)) - 1;

      return (s >= lb) && (s <= ub);
    }

    // Checks if the value of an usigned integer fits into the lowest 'B' bits:
    template <unsigned B, std::unsigned_integral U>
    bool check(U u)
    {
      static_assert(B > 0);
      static_assert(B <= 8 * sizeof(U));

      static constexpr U ub = (static_cast<U>(1) << (B - 1)) - 1;

      return (u <= ub);
    }

    // Converts an integer to 'U' while preserving its bit pattern, masks 'B'
    // lowest bytes, and left-shifts the result by 'S':
    template <unsigned_integral U, unsigned B, unsigned S = 0, std::integral I>
    constexpr U convert_mask_shift(I i)
    {
      assert(check<B>(i));

      static_assert(B > 0);
      static_assert(B <= 8 * sizeof(U));
      static_assert(B <= 8 * sizeof(I));

      using UI = std::make_unsigned_t<I>; // no effect if 'I' is unsigned

      // The following conversion preserves bit pattern since C++20, or, if
      // signed integral types uses two's complement representation:
      UI ui = static_cast<UI>(i); 

      // Check bit pattern preservation:
      static_assert(sizeof(I) == sizeof(UI));
      assert(std::memcmp(&i, &ui, sizeof(I)) == 0);

      // Mask B lowest bits:
      ui &= mask<UI, B>();

      // Convert to the resulting type:
      U u = static_cast<U>(ui);

      // Back-check the reverse conversion:
      assert(static_cast<UI>(u) == ui);

      static_assert(B + S <= 8 * sizeof(U));

      // Shift needs to be a multiple of 'B':
      static_assert(S % B == 0);

      return u << S;
    }

    template <unsigned B, unsigned S, unsigned_integral U, std::integral Head,
              std::integral... Tail>
    U pack_ints_helper(Head head, Tail... tail)
    {
      U part = intpack::convert_mask_shift<U, B, S>(head);

      if constexpr (sizeof...(Tail) == 0)
        return part;
      else
        return part | pack_ints_helper<B, S + B, U>(tail...);
    }

  } // namespace intpack

  template <intpack::unsigned_integral U, unsigned B = 0, std::integral... Ints>
  U pack_ints(Ints... ints)
  {
    static constexpr unsigned BB =
        (B > 0) ? B : sizeof(U) * 8 / (sizeof...(Ints));

    static_assert(BB > 0);

    return intpack::pack_ints_helper<BB, 0, U>(ints...);
  }

} // namespace utils

#endif
