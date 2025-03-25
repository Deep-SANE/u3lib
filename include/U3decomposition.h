#ifndef U3LIB_DECOMPOSITION_H
#define U3LIB_DECOMPOSITION_H

#include <algorithm>
#include <utility>
#include <vector>

#include "U3basis.h"

namespace u3
{

  std::pair<std::vector<u3::irrep>, std::vector<int>>
  decomposition(const u3::irrep& u31, const u3::irrep& u32)
  {
    std::vector<u3::irrep> u3cs;
    std::vector<int> rhomaxs;

    u3::irrep u3a, u3b;

    if (u31.dim() <= u32.dim())
    {
      u3a = u31;
      u3b = u32;
    }
    else
    {
      u3a = u32;
      u3b = u31;
    }

    for (int n12 = u3a.n13(); n12 >= u3a.n23(); n12--)
    {
      for (int n22 = u3a.n23(); n22 >= u3a.n33(); n22--)
      {
        for (int n11 = n12; n11 >= n22; n11--)
        {
          int b13 = u3a.n13() - n12, b23 = u3a.n23() - n22, b33 = u3a.n33(),
              b12 = n12 - n11, b22 = n22, b11 = n11;

          u3::irrep u3c = u3b;

          u3c.n33(u3c.n33() + b13);
          if (!u3c.check())
            continue;

          u3c.n23(u3c.n23() + b12);
          if (!u3c.check())
            continue;

          u3c.n13(u3c.n13() + b11);
          if (!u3c.check())
            continue;

          u3c.n33(u3c.n33() + b23);
          if (!u3c.check())
            continue;

          u3c.n23(u3c.n23() + b22);
          if (!u3c.check())
            continue;

          u3c.n33(u3c.n33() + b33);
          if (!u3c.check())
            continue;

          auto iter = std::find(u3cs.begin(), u3cs.end(), u3c);
          if (iter == u3cs.end())
          {
            u3cs.push_back(u3c);
            rhomaxs.push_back(1);
          }
          else
          {
            auto index = std::distance(u3cs.begin(), iter);
            rhomaxs[index]++;
          }
        }
      }
    }

    return std::make_pair(std::move(u3cs), std::move(rhomaxs));
  }

  // compute the outer multiplicity directly, before normalization
  int Rho(const u3::irrep& u3a, const u3::irrep& u3b, const u3::irrep& u3c)
  {
    int rho = 0, n11 = u3c.n13() - u3b.n13();

    for (int n12 = u3a.n13(); n12 >= std::max(u3a.n23(), n11); n12--)
    {
      for (int n22 = std::min(u3a.n23(), n11); n22 >= u3a.n33(); n22--)
      {
        int b13 = u3a.n13() - n12, b23 = u3a.n23() - n22, b33 = u3a.n33(),
            b12 = n12 - n11, b22 = n22, b11 = n11;

        u3::irrep u3C = u3b;

        u3C.n33(u3C.n33() + b13);
        if (!u3C.check())
          continue;

        u3C.n23(u3C.n23() + b12);
        if (!u3C.check())
          continue;

        u3C.n13(u3C.n13() + b11);
        if (!u3C.check())
          continue;

        u3C.n33(u3C.n33() + b23);
        if (!u3C.check())
          continue;

        u3C.n23(u3C.n23() + b22);
        if (!u3C.check())
          continue;

        u3C.n33(u3C.n33() + b33);
        if (!u3C.check())
          continue;

        if (u3C == u3c)
          rho++;
      }
    }

    return rho;
  }

} // namespace u3

#endif
