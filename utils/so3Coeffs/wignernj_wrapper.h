#ifndef UTILS_WIGNERNJ_WRAPPER_H
#define UTILS_WIGNERNJ_WRAPPER_H

#include <cassert>

#ifdef UTILS_WIGNERNJ_LOG
#include <fstream>
#endif

// libwignernj C++ header (header-only inline wrappers over the C library).
// Defines the global namespace ::wignernj; we reach it with a leading "::" to
// avoid clashing with the utils::wignernj class defined below.
#include <wignernj.hpp>

namespace utils
{
  // Drop-in replacement for utils::wigxjpf, exposing the same static
  // interface (wig3jj / wig6jj / wig9jj + lifecycle hooks) so that
  // utils::wigner can alias either backend.
  //
  // Unlike wigxjpf, libwignernj needs no global symbol table and no explicit
  // per-thread temporaries: every per-thread cache is lazily allocated in
  // thread-local storage and is therefore safe under OpenMP without setup.
  // The one shared piece of state is the prime table (primes_init), which is
  // built on the first symbol evaluation and is NOT thread-safe; initialize()
  // forces it once while still serial, before any parallel region.
  class wignernj
  {
  public:
    // Call once in the main thread before spawning any workers. A trivial
    // 3j evaluation (which equals 1) triggers the one-time, non-thread-safe
    // prime-table initialization deterministically while serial.
    static void initialize()
    {
      (void) ::wignernj::symbol3j<double>(0, 0, 0, 0, 0, 0);
    }

    // No per-thread setup is required: caches are lazy and thread-local.
    static void initialize_thread() {}

    // No per-thread teardown is required. Caches are reclaimed when the
    // thread exits (or via finalize()); releasing them at the end of every
    // parallel region would only force re-allocation in the next one.
    static void finalize_thread() {}

    // Release the calling (main) thread's caches at shutdown.
    static void finalize()
    {
      ::wignernj::thread_cleanup();
    }

    // 3j / 6j / 9j symbols. Arguments are 2*j (and 2*m) integers, matching the
    // wigxjpf convention exactly (verified numerically against wigxjpf over
    // 6e5 randomized triples: identical to machine precision).
    static auto wig3jj(int a, int b, int c, int d, int e, int f)
    {
#ifdef UTILS_WIGNERNJ_LOG
      static std::ofstream f_log("wignernj_wig3jj.log");
      f_log << a << " " << b << " " << c << " " << d << " " << e << " " << f
            << "\n";
#endif
      return ::wignernj::symbol3j<double>(a, b, c, d, e, f);
    }

    static auto wig6jj(int a, int b, int c, int d, int e, int f)
    {
#ifdef UTILS_WIGNERNJ_LOG
      static std::ofstream f_log("wignernj_wig6jj.log");
      f_log << a << " " << b << " " << c << " " << d << " " << e << " " << f
            << "\n";
#endif
      return ::wignernj::symbol6j<double>(a, b, c, d, e, f);
    }

    static auto wig9jj(int a, int b, int c, int d, int e, int f, int g, int h,
                       int i)
    {
#ifdef UTILS_WIGNERNJ_LOG
      static std::ofstream f_log("wignernj_wig9jj.log");
      f_log << a << " " << b << " " << c << " " << d << " " << e << " " << f
            << " " << g << " " << h << " " << i << "\n";
#endif
      return ::wignernj::symbol9j<double>(a, b, c, d, e, f, g, h, i);
    }
  };
} // namespace utils

#endif
