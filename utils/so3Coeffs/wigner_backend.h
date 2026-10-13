#ifndef UTILS_WIGNER_BACKEND_H
#define UTILS_WIGNER_BACKEND_H

// Exposes the Wigner-symbol backend under the backend-neutral name
// utils::wigner. libwignernj is the sole backend; the selecting macro
// USE_WIGNERNJ is emitted by config/make.d/wignernj.mk.

#if defined(USE_WIGNERNJ)
#include <so3Coeffs/wignernj_wrapper.h>
#else
#error "No Wigner backend selected: compile with -DUSE_WIGNERNJ"
#endif

namespace utils
{
  using wigner = wignernj;
} // namespace utils

#endif
