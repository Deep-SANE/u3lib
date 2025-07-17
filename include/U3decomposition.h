#ifndef U3LIB_DECOMPOSITION_H
#define U3LIB_DECOMPOSITION_H

#include <algorithm>
#include <utility>
#include <vector>

#include "U3irrep.h"

namespace su3
{
  // original Tomas' version of 'su3mult' (renamed to 'Rho' here):
  inline int Rho(int lm1, int mu1, int lm2, int mu2, int lm3, int mu3)
  {
    int result = 0;
    int phase = lm1 + lm2 - lm3 - mu1 - mu2 + mu3;

    if (phase % 3 == 0)
    {
      int l1, l2, l3, m1, m2, m3;
      int redphase = phase / 3;

      if (phase >= 0)
      {
        l1 = lm1;
        l2 = lm2;
        l3 = lm3;
        m1 = mu1;
        m2 = mu2;
        m3 = mu3;
      }
      else
      {
        l1 = mu1;
        l2 = mu2;
        l3 = mu3;
        m1 = lm1;
        m2 = lm2;
        m3 = lm3;

        redphase = -redphase;
      }

      phase = redphase + m1 + m2 - m3;

      int invl1 = std::min(l1 - redphase, m2);
      if (invl1 < 0)
        return result;

      int invl2 = std::min(l2 - redphase, m1);
      if (invl2 < 0)
        return result;

      result =
          std::max(std::min(phase, invl2) - std::max(phase - invl1, 0) + 1, 0);
    }

    return result;
  }

  // version for 'su3::irrep's:
  inline int Rho(const su3::irrep& su3a, const su3::irrep& su3b,
                 const su3::irrep& su3c)
  {
    return su3::Rho(su3a.lam(), su3a.mu(), su3b.lam(), su3b.mu(), su3c.lam(),
                    su3c.mu());
  }

  inline std::pair<std::vector<su3::irrep>, std::vector<int>>
  decomposition(const su3::irrep& su3a, const su3::irrep& su3b)
  {
    std::vector<su3::irrep> su3cs;
    std::vector<int> rhomaxs;

    auto lm1 = su3a.lam();
    auto mu1 = su3a.mu();

    auto lm2 = su3b.lam();
    auto mu2 = su3b.mu();

    su3cs.reserve(3 * (lm1 + mu1 + lm2 + mu2));
    rhomaxs.reserve(3 * (lm1 + mu1 + lm2 + mu2));

    int ix_max = std::min(mu1, lm2);
    for (int ix = 0; ix <= ix_max; ++ix)
    {
      int iy_max = std::min(lm1, lm2 + mu2 - ix);
      for (int iy = 0; iy <= iy_max; ++iy)
      {
        int iz_min = std::max(0, mu2 - mu1 + ix - iy);
        int iz_max = std::min(mu2, lm2 + mu2 - ix - iy);
        for (int iz = iz_min; iz <= iz_max; ++iz)
        {
          int lm3 = lm1 + lm2 + mu2 - 2 * iy - ix - iz;
          int mu3 = mu1 - mu2 - ix + iy + 2 * iz;
          int ip = (2 * (mu3 - mu1) + lm3 - lm1 - lm2 + mu2) / 3;
          int in = (2 * (lm1 + lm2 - lm3) + mu1 + mu2 - mu3) / 3;
          int irho = ix - std::max({0, -ip, in - lm1}) + 1;

          if (irho != 1)
            continue;

          rhomaxs.push_back(su3::Rho(lm1, mu1, lm2, mu2, lm3, mu3));
          su3cs.emplace_back(lm3, mu3);
        }
      }
    }

    return std::make_pair(std::move(su3cs), std::move(rhomaxs));
  }
} // namespace su3

namespace u3
{
  inline int Rho(const u3::irrep& u3a, const u3::irrep& u3b,
                 const u3::irrep& u3c)
  {
    return su3::Rho(u3a, u3b, u3c);
  }

  inline std::pair<std::vector<u3::irrep>, std::vector<int>>
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

          u3c.normalize();
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
} // namespace u3

#endif
