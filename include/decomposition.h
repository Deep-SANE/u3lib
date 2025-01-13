// protection against multiple header file inclusion
#ifndef DECOMPOSITION_H
#define DECOMPOSITION_H

#include <cmath>
#include <iostream>
#include <vector>

#include "U3basis.h"

namespace u3
{

  std::vector<u3::irrep> decomposition(u3::irrep& u31, u3::irrep& u32,
                                       std::vector<int>& rhomaxs)
  {

    std::vector<u3::irrep> u3cs;

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

    for (int n12 = u3a.n13; n12 >= u3a.n23; n12--)
    {
      for (int n22 = u3a.n23; n22 >= u3a.n33; n22--)
      {
        for (int n11 = n12; n11 >= n22; n11--)
        {

          int b13 = u3a.n13 - n12, b23 = u3a.n23 - n22, b33 = u3a.n33,
              b12 = n12 - n11, b22 = n22, b11 = n11;

          u3::irrep u3c = u3b;

          u3c.n33 += b13;
          if (!u3c.check())
            continue;

          u3c.n23 += b12;
          if (!u3c.check())
            continue;

          u3c.n13 += b11;
          if (!u3c.check())
            continue;

          u3c.n33 += b23;
          if (!u3c.check())
            continue;

          u3c.n23 += b22;
          if (!u3c.check())
            continue;

          u3c.n33 += b33;
          if (!u3c.check())
            continue;

          if (u3cs.size() == 0)
          {
            u3cs.push_back(u3c);
            rhomaxs.push_back(1);
          }
          else
          {
            int flag = 0;
            for (int i = 0; i < u3cs.size(); i++)
            {
              if (u3c.n13 == u3cs[i].n13 && u3c.n23 == u3cs[i].n23 &&
                  u3c.n33 == u3cs[i].n33)
              {
                rhomaxs[i]++;
                flag = 1;
              }
            }

            if (flag == 0)
            {
              u3cs.push_back(u3c);
              rhomaxs.push_back(1);
            }
          }
        }
      }
    }

    return u3cs;
  }

} // namespace u3

#endif