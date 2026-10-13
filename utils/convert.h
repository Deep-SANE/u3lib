#ifndef UTILS_CONVERT_H
#define UTILS_CONVERT_H

#include <cassert>
#include <cstdint>
#include <limits>

namespace utils
{
  // helper function for safe-asserted conversion from `int` to `uint8_t`:
  inline uint8_t itou8(int i)
  {
    static const int max =
        static_cast<int>(std::numeric_limits<uint8_t>::max());

    assert(i >= 0);
    assert(i <= max);

    return static_cast<uint8_t>(i);
  }

  // helper function for safe-asserted conversion from `int` to `int8_t`:
  inline int8_t ito8(int i)
  {
    static const int max =
        static_cast<int>(std::numeric_limits<int8_t>::max());
    static const int min =
        static_cast<int>(std::numeric_limits<int8_t>::min());

    assert(i >= min);
    assert(i <= max);

    return static_cast<int8_t>(i);
  }

  // helper function for safe-asserted conversion from `int` to `uint16_t`:
  inline uint16_t itou16(int i)
  {
    static const int max =
        static_cast<int>(std::numeric_limits<uint16_t>::max());

    assert(i >= 0);
    assert(i <= max);

    return static_cast<uint16_t>(i);
  }

  // helper function for safe-asserted conversion from `int` to `int16_t`:
  inline int16_t ito16(int i)
  {
    static const int max =
        static_cast<int>(std::numeric_limits<int16_t>::max());
    static const int min =
        static_cast<int>(std::numeric_limits<int16_t>::min());

    assert(i >= min);
    assert(i <= max);

    return static_cast<int16_t>(i);
  }

} // namespace utils

#endif
