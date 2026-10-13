#ifndef UTILS_PRINT_H
#define UTILS_PRINT_H

/*
// C++23:
#if defined __cpp_lib_print && defined __cpp_lib_format &&                     \
    defined __cpp_lib_format_ranges

#include <print>

using std::print;
// using std::println;
*/

// {fmt} library:
//#elif defined(HAVE_FMT)
#if defined(HAVE_FMT)

#include <fmt/format.h>
#include <fmt/ranges.h>

namespace utils
{
  using fmt::print;
}

// none of the above options:
#else

#include <iostream>

namespace utils
{
  template <typename... Ts>
  void print(const Ts&...)
  {
    std::cout << "{fmt} library not configured; this is a fallback message!\n";
  }
} // namespace utils

#endif

#endif
