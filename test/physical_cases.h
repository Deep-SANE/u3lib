#pragma once
#include "test_helpers.h"
namespace test
{
  inline void physical(Suite& suite, Reference& ref)
  {
    for (const auto& x : couplings())
      suite.run(
          "Physical coefficients " + label(x) +
              " against documented values and orthogonality",
          [&]
          {
            auto pa = u3::physicalBasis(x.a), pb = u3::physicalBasis(x.b),
                 pc = u3::physicalBasis(x.c);
            int rho = u3::Rho(x.a, x.b, x.c);
            auto [indices, coeffs] = u3::WigU3SO3(rho, pa, pb, pc);
            require(!coeffs.empty(), "Empty physical coefficient block");
            // Independent contraction for one complete kappa/rho block at
            // M_C=0. This exercises a different magnetic projection from the
            // library's highest-M reduced-coefficient construction.
            {
              auto [ai, bi, ci] = indices[indices.size() / 2];
              int LA = pa[ai].L(), LB = pb[bi].L(), LC = pc[ci].L(), MA = 0;
              double angular = 0;
              for (int m = -std::min(LA, LB); m <= std::min(LA, LB); ++m)
              {
                double candidate = qcx_utils::cgcs().get(
                    {2 * LA, 2 * LB, 2 * LC, 2 * m, -2 * m, 0});
                if (std::abs(candidate) > std::abs(angular))
                {
                  angular = candidate;
                  MA = m;
                }
              }
              require(std::abs(angular) > 1e-12,
                      "Nonzero angular CG for contraction");
              auto ga = u3::gelfandBasis(x.a), gb = u3::gelfandBasis(x.b),
                   gc = u3::gelfandBasis(x.c);
              auto [ia, ta] = u3::phys_to_cano(pa[ai], MA, ga);
              auto [ib, tb] = u3::phys_to_cano(pb[bi], -MA, gb);
              auto [ic, tc] = u3::phys_to_cano(pc[ci], 0, gc);
              auto canonical = u3::WigU3U2(x.a, x.b, x.c, rho, -1, ga, gb, gc);
              auto actual = coeffs[indices.size() / 2];
              std::vector<double> contracted(actual.size(), 0);
              // Explicit product-basis sum, with full canonical SU(2) CG
              // factors.
              for (int c = 0; c < int(ic.size()); ++c)
              {
                auto g = gc[ic[c]];
                auto allowed = u3::expand_to_uncoupled(g, ga, gb);
                for (auto [aidx, bidx] : allowed)
                {
                  auto ita = std::find(ia.begin(), ia.end(), aidx),
                       itb = std::find(ib.begin(), ib.end(), bidx);
                  if (ita == ia.end() || itb == ib.end())
                    continue;
                  int a = std::distance(ia.begin(), ita),
                      b = std::distance(ib.begin(), itb);
                  auto ag = ga[aidx], bg = gb[bidx];
                  auto acc = canonical.find_row(ag.n12(), ag.n22(), bg.n12(),
                                                bg.n22(), g.n12(), g.n22());
                  if (!acc.allowed)
                    continue;
                  double cg = qcx_utils::cgcs().get(
                      {ag.Lambdax2(), bg.Lambdax2(), g.Lambdax2(),
                       ag.MLambdax2(), bg.MLambdax2(), g.MLambdax2()});
                  for (int ka = 0; ka < ta.cols(); ++ka)
                    for (int kb = 0; kb < tb.cols(); ++kb)
                      for (int kc = 0; kc < tc.cols(); ++kc)
                        for (int r = 0; r < rho; ++r)
                        {
                          int index =
                              ((ka * tb.cols() + kb) * tc.cols() + kc) * rho +
                              r;
                          contracted[index] += ta(a, ka) * tb(b, kb) *
                                               tc(c, kc) * cg *
                                               canonical.get(acc, r);
                        }
                }
              }
              for (int i = 0; i < int(actual.size()); ++i)
                near(contracted[i], angular * actual[i],
                     "Independent physical/canonical contraction", 1e-10,
                     1e-10);
            }
            for (int z = 0; z < int(pc.size()); ++z)
            {
              int nk = pc[z].kappa();
              Eigen::MatrixXd gram = Eigen::MatrixXd::Zero(nk * rho, nk * rho);
              for (int t = 0; t < int(indices.size()); ++t)
              {
                auto [a, b, c] = indices[t];
                if (c != z)
                  continue;
                const auto& v = coeffs[t];
                for (int ka = 0; ka < pa[a].kappa(); ++ka)
                  for (int kb = 0; kb < pb[b].kappa(); ++kb)
                  {
                    Eigen::VectorXd row(nk * rho);
                    for (int kc = 0; kc < nk; ++kc)
                      for (int r = 0; r < rho; ++r)
                      {
                        double value = v(ka, kb, kc, r);
                        row[kc * rho + r] = value;
                        ref.check(key(label(x), pa[a].L(), pb[b].L(), pc[c].L(),
                                      ka, kb, kc, r),
                                  value);
                        if (x.a.lam() == 0 && x.a.mu() == 0)
                          near(value, kb == kc ? 1 : 0,
                               "Analytic scalar physical coefficient");
                      }
                    gram += row * row.transpose();
                  }
              }
              near(
                  (gram - Eigen::MatrixXd::Identity(nk * rho, nk * rho)).norm(),
                  0, "Physical reduced-CG orthogonality", 1e-9, 0);
            }
          });
  }
} // namespace test
