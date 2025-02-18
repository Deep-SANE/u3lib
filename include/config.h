#ifndef U3LIB_CONFIG_H
#define U3LIB_CONFIG_H

#include "gelfand/gelfand.h"
#include "irrep/irrep.h"

namespace u3
{
  using irrep = internal::irrep_int_rep_t;
  using gelfand = internal::gelfand_int_rep_t;
} // namespace u3

#endif
