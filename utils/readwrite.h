#ifndef UTILS_READWRITE_H
#define UTILS_READWRITE_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

// readwrite.h itself no longer needs Eigen: the sparse-matrix and
// eigen-solution file I/O moved to SAM/utils/readwrite.h. the includes stay
// because several u3lib/u4lib/SAM files use Eigen without including it
// directly and rely on getting it through this header.
#include <Eigen/Dense>
#include <Eigen/Sparse>

// largest number of bytes passed to one stream write/read call. macOS fails
// any single write(2)/read(2) above 2^31 - 1 bytes (e.g. the 3.3 GB values of
// a CSR matrix), so large vectors go through in pieces of this size; the
// bytes in the file are unchanged.
#ifndef UTILS_RW_CHUNK_BYTES
#define UTILS_RW_CHUNK_BYTES (std::size_t(1) << 30)
#endif

namespace utils
{
  // helper function for binary writing a single value to output stream:
  template <typename T>
  inline void writeBinary(std::ostream& s, const T& t)
  {
    const char* ptr = reinterpret_cast<const char*>(&t);
    s.write(ptr, sizeof(t));
  }

  // helper function for binary writing a vector to output stream:
  template <typename T>
  inline void writeBinary(std::ostream& s, const std::vector<T>& v)
  {
    const char* ptr = reinterpret_cast<const char*>(v.data());
    std::size_t left = v.size() * sizeof(T);
    while (left > 0 && s)
    {
      const std::size_t n = std::min(left, std::size_t(UTILS_RW_CHUNK_BYTES));
      s.write(ptr, static_cast<std::streamsize>(n));
      ptr += n;
      left -= n;
    }
  }

  // helper function for binary reading a single value from input stream:
  template <typename T>
  inline void readBinary(std::istream& s, T& t)
  {
    char* ptr = reinterpret_cast<char*>(&t);
    s.read(ptr, sizeof(t));
  }

  // helper function for binary reading a vector from input stream:
  template <typename T>
  inline void readBinary(std::istream& s, std::vector<T>& v)
  {
    assert(v.size() > 0);
    char* ptr = reinterpret_cast<char*>(v.data());
    std::size_t left = v.size() * sizeof(T);
    while (left > 0 && s)
    {
      const std::size_t n = std::min(left, std::size_t(UTILS_RW_CHUNK_BYTES));
      s.read(ptr, static_cast<std::streamsize>(n));
      ptr += n;
      left -= n;
    }
  }

  /*
  // helper function for binary writing a single value to output stream:
  template <typename T, typename V>
  inline void write(std::ostream& s, const V& v)
  {
    T t = v; // 'T' is the written type
    const char* ptr = reinterpret_cast<const char*>(&t);
    s.write(ptr, sizeof(t));
  }

  // helper function for binary reading a single value from input stream:
  template <typename T, typename V>
  inline void read(std::istream& s, V& v)
  {
    T t; // 'T' is the read type
    char* ptr = reinterpret_cast<char*>(&t);
    s.read(ptr, sizeof(t));
    v = t;
  }
  */
} // namespace utils

#endif
