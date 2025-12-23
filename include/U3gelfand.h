#ifndef U3LIB_GELFAND_H
#define U3LIB_GELFAND_H

#include <vector>

#include "internal/gelfand/gelfand.h"

namespace u3
{
  using gelfand = internal::gelfand_int_rep_t;
  using gelfands = std::vector<u3::gelfand>;
}

#endif
