#ifndef UTILS_H
#define UTILS_H

#include <so3Coeffs/wigner_backend.h>
#include <so3Coeffs/SO3_CGC.h>
#include <so3Coeffs/SO3_Racah.h>
#include <so3Coeffs/SO3_Fano.h>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace utils
{
  // Process-global init/finalize are idempotent: with the u3 -> utils and
  // u4 -> u3 -> utils cascades, initialize() may be reached more than once,
  // but MPI_Init (and the SO3 caches) must run exactly once.
  // The first caller supplies argc/argv; later calls are no-ops.
  inline bool& initialized_flag()
  {
    static bool flag = false;
    return flag;
  }

  inline void initialize(int* argc, char*** argv)
  {
    if (initialized_flag()) return;
    initialized_flag() = true;

    wigner::initialize();

#ifdef UTILS_SO3_CGS_CACHE
    SO3_CGC::initialize_cache();
#endif
#ifdef UTILS_SO3_RACAH_CACHE
    SO3_Racah::initialize_cache();
#endif
#ifdef UTILS_SO3_FANO_CACHE
    SO3_Fano::initialize_cache();
#endif

#ifdef USE_MPI
    // THREAD_SERIALIZED: the process is multithreaded (OpenMP) and MPI may be
    // called from a worker thread, but never by two threads at once. The expand
    // dynamic work-queue relies on this (it claims tensor batches with an RMA
    // atomic from inside an `omp single`). A superset of the old plain MPI_Init
    // (== THREAD_SINGLE): every existing call site is outside parallel regions,
    // so it stays valid. Modern MPIs always provide at least SERIALIZED.
    int mpi_thread_provided = 0;
    MPI_Init_thread(argc, argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
#endif
  }

  // Per-thread init/finalize are idempotent per thread (the u3/u4 _thread
  // cascades plus legacy explicit calls can reach them more than once on the
  // same thread).
  inline bool& thread_initialized_flag()
  {
    static thread_local bool flag = false;
    return flag;
  }

  inline void initialize_thread()
  {
    if (thread_initialized_flag()) return;
    thread_initialized_flag() = true;

    wigner::initialize_thread();

#ifdef USE_OPENMP
    if (omp_get_thread_num() != 0)
    {
#ifdef UTILS_SO3_CGS_CACHE
      SO3_CGC::initialize_cache();
#endif
#ifdef UTILS_SO3_RACAH_CACHE
      SO3_Racah::initialize_cache();
#endif
#ifdef UTILS_SO3_FANO_CACHE
      SO3_Fano::initialize_cache();
#endif
    }
#endif
  }

  inline void finalize()
  {
    if (!initialized_flag()) return;
    initialized_flag() = false;

#ifdef USE_MPI
    MPI_Finalize();
#endif

#ifdef UTILS_SO3_FANO_CACHE
    SO3_Fano::finalize_cache();
#endif
#ifdef UTILS_SO3_RACAH_CACHE
    SO3_Racah::finalize_cache();
#endif
#ifdef UTILS_SO3_CGS_CACHE
    SO3_CGC::finalize_cache();
#endif

    wigner::finalize();
  }

  inline void finalize_thread()
  {
    if (!thread_initialized_flag()) return;
    thread_initialized_flag() = false;

#ifdef USE_OPENMP
    if (omp_get_thread_num() != 0)
    {
#ifdef UTILS_SO3_FANO_CACHE
      SO3_Fano::finalize_cache();
#endif
#ifdef UTILS_SO3_RACAH_CACHE
      SO3_Racah::finalize_cache();
#endif
#ifdef UTILS_SO3_CGS_CACHE
      SO3_CGC::finalize_cache();
#endif
    }
#endif

    wigner::finalize_thread();
  }

} // namespace utils

#endif
