#ifndef U3LIB_U3LIB_H
#define U3LIB_U3LIB_H

#include <utils.h>   // foundation (Wigner backend + SO3 caches); folded in below

#include <U3U2xU3U2.h>
#include <U3SO3xU3SO3.h>
#include <U3_U6Z6Coeffs.h>
#include <U3_U9Coeffs.h>

namespace u3
{
  // initialize()/finalize() also set up the utils foundation, so callers
  // (u4lib, SAM, or standalone users) only need to init the highest layer they
  // use. argc/argv are forwarded to utils::initialize (MPI); they
  // default to null for non-MPI / standalone use.
  //
  // All four are idempotent: because u4 cascades u4 -> u3 -> utils, and
  // legacy programs still call utils::initialize()/u3::initialize()
  // explicitly, u3::initialize() may be reached more than once. The guards make
  // repeated calls no-ops (process-wide for init/finalize, per-thread for the
  // _thread variants).
  inline bool& initialized_flag() { static bool f = false; return f; }
  inline bool& thread_initialized_flag() { static thread_local bool f = false; return f; }

  // Read-only query of the process-global initialization state.
  inline bool is_initialized() { return initialized_flag(); }

  inline void initialize(int* argc = nullptr, char*** argv = nullptr)
  {
    if (initialized_flag()) return;
    initialized_flag() = true;
    utils::initialize(argc, argv);
    WigU3U2_db::initialize();
    WigU3SO3_db::initialize();
    U3U6_db::initialize();
    U3Z6_db::initialize();
    U3U9_db::initialize();
  }

  inline void finalize()
  {
    if (!initialized_flag()) return;
    initialized_flag() = false;
    WigU3U2_db::finalize();
    WigU3SO3_db::finalize();
    U3U6_db::finalize();
    U3Z6_db::finalize();
    U3U9_db::finalize();
    utils::finalize();   // foundation torn down last
  }

  inline void initialize_thread()
  {
    if (thread_initialized_flag()) return;
    thread_initialized_flag() = true;
    utils::initialize_thread();
    WigU3U2_db::initialize_thread();
    WigU3SO3_db::initialize_thread();
    U3U6_db::initialize_thread();
    U3Z6_db::initialize_thread();
    U3U9_db::initialize_thread();
  }

  inline void finalize_thread()
  {
    if (!thread_initialized_flag()) return;
    thread_initialized_flag() = false;
    WigU3U2_db::finalize_thread();
    WigU3SO3_db::finalize_thread();
    U3U6_db::finalize_thread();
    U3Z6_db::finalize_thread();
    U3U9_db::finalize_thread();
    utils::finalize_thread();
  }

  inline void clear()
  {
    wigU3U2db().clear();
    wigU3SO3db().clear();
    u3u6db().clear();
    u3z6db().clear();
    u3u9db().clear();
  }
}

#endif
