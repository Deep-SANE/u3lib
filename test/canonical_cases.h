#pragma once
#include "test_helpers.h"
namespace test
{
  inline void canonical(Suite& suite, Reference& ref)
  {
    for (const auto& x : couplings())
      suite.run(
          "Canonical coefficients " + label(x) +
              " against documented values and normalization",
          [&]
          {
            auto ga = u3::gelfandBasis(x.a), gb = u3::gelfandBasis(x.b),
                 gc = u3::gelfandBasis(x.c);
            int rho = u3::Rho(x.a, x.b, x.c);
            require(rho > 0, "Forbidden test coupling");
            auto coeff = u3::WigU3U2(x.a, x.b, x.c, rho, -1, ga, gb, gc);
            // One representative per U(2) irrep; store every allowed reduced
            // row.
            for (auto a : ga)
              if (a.n11() == a.n12())
                for (auto b : gb)
                  if (b.n11() == b.n12())
                    for (auto c : gc)
                      if (c.n11() == c.n12())
                      {
                        auto acc = coeff.find_row(a.n12(), a.n22(), b.n12(),
                                                  b.n22(), c.n12(), c.n22());
                        if (!acc.allowed)
                          continue;
                        for (int r = 0; r < rho; ++r)
                          ref.check(key(label(x), int(a.n12()), int(a.n22()),
                                        int(b.n12()), int(b.n22()),
                                        int(c.n12()), int(c.n22()), r),
                                    coeff.get(acc, r));
                      }
            // Reconstruct full canonical CG columns for EVERY output Gelfand
            // state.
            for (auto c : gc)
            {
              auto ids = u3::expand_to_uncoupled(c, ga, gb);
              Eigen::MatrixXd v = Eigen::MatrixXd::Zero(ids.size(), rho);
              for (int i = 0; i < int(ids.size()); ++i)
              {
                auto a = ga[ids[i].first], b = gb[ids[i].second];
                auto acc = coeff.find_row(a.n12(), a.n22(), b.n12(), b.n22(),
                                          c.n12(), c.n22());
                if (!acc.allowed)
                  continue;
                double cg = qcx_utils::cgcs().get(
                    {a.Lambdax2(), b.Lambdax2(), c.Lambdax2(), a.MLambdax2(),
                     b.MLambdax2(), c.MLambdax2()});
                v.row(i) = cg * coeff.get_row(acc);
                if (x.c.lam() == 0 && x.c.mu() == 0)
                  near(v(i, 0),
                       a.conjugate() == b
                           ? a.conjugate_phase() / std::sqrt(double(ga.size()))
                           : 0,
                       "Analytic singlet phase");
                if (x.a.lam() == 0 && x.a.mu() == 0)
                  near(v(i, 0), 1, "Analytic scalar identity");
              }
              orthogonal(v, "Full canonical CG orthonormality");
            }
          });
  }
} // namespace test
