#ifndef UTILS_COUNTER_H
#define UTILS_COUNTER_H

#include <cstddef>
#include <source_location>

#include <print.h>

namespace utils
{
  class Counter
  {
    std::source_location srcl_;
    std::size_t cnt_;

  public:
    Counter(const std::source_location srcl = std::source_location::current())
        : srcl_(srcl), cnt_(0)
    { }

    ~Counter()
    {
      print("{}: {}: {}\n", srcl_.file_name(), srcl_.function_name(), cnt_);
    }

    Counter& operator++()
    {
      ++cnt_; 
      return *this;
    }
  };

} // namespace utils

#endif
