#ifndef U3LIB_U3LIB_H
#define U3LIB_U3LIB_H

#include <U3U2xU3U2.h>
#include <U3SO3xU3SO3.h>
#include <U3_U6Z6Coeffs.h>
#include <U3_U9Coeffs.h>

namespace u3
{
  inline void initialize()
  {
    WigU3U2_db::initialize();
    WigU3SO3_db::initialize();
    U3U6_db::initialize();
    U3Z6_db::initialize();
    U3U9_db::initialize();
  }

  inline void finalize()
  {
    WigU3U2_db::finalize();
    WigU3SO3_db::finalize();
    U3U6_db::finalize();
    U3Z6_db::finalize();
    U3U9_db::finalize();
  }

  inline void initialize_thread()
  {
    WigU3U2_db::initialize_thread();
    WigU3SO3_db::initialize_thread();
    U3U6_db::initialize_thread();
    U3Z6_db::initialize_thread();
    U3U9_db::initialize_thread();
  }

  inline void finalize_thread()
  {
    WigU3U2_db::finalize_thread();
    WigU3SO3_db::finalize_thread();
    U3U6_db::finalize_thread();
    U3Z6_db::finalize_thread();
    U3U9_db::finalize_thread();
  }
}

#endif
