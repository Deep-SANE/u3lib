#ifndef U3LIB_U3LIB_H
#define U3LIB_U3LIB_H

#include <U3U2xU3U2.h>

namespace u3
{

  inline void initialize()
  {
    WigU3U2_db::initialize();
  }

  inline void finalize()
  {
    WigU3U2_db::finalize();
  }

#ifdef USE_OPENMP
  inline void initialize_thread()
  {
    WigU3U2_db::initialize_thread();
  }

  inline void finalize_thread()
  {
    WigU3U2_db::finalize_thread();
  }
#endif
}

#endif
