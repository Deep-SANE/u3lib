#pragma once
#include <U3generatorsRME.h>
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
                double candidate = utils::cgcs().get(
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
                  double cg = utils::cgcs().get(
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

    suite.run(
        "Physical generator RMEs (L, Q) zero-based multiplicity reporting",
        [&]
        {
          u3::physicalvector p22_2(u3::irrep(2, 2), 2);
          require(p22_2.kappa() == 2, "kappa(2,2; L=2) == 2");
          std::vector<int> bk_L, kk_L;
          auto rme_L = u3::RME_L(p22_2, p22_2, bk_L, kk_L);
          require(bk_L.size() == 2 && kk_L.size() == 2, "RME_L out-param size");
          require(bk_L[0] == 0 && bk_L[1] == 1, "RME_L bra_kappa zero-based");
          require(kk_L[0] == 0 && kk_L[1] == 1, "RME_L ket_kappa zero-based");
          require(rme_L.rows() == 2 && rme_L.cols() == 2, "RME_L shape");

          std::vector<int> bk_Q, kk_Q;
          auto rme_Q = u3::RME_Q(p22_2, p22_2, bk_Q, kk_Q);
          require(bk_Q.size() == 2 && kk_Q.size() == 2, "RME_Q out-param size");
          require(bk_Q[0] == 0 && bk_Q[1] == 1, "RME_Q bra_kappa zero-based");
          require(kk_Q[0] == 0 && kk_Q[1] == 1, "RME_Q ket_kappa zero-based");
          require(rme_Q.rows() == 2 && rme_Q.cols() == 2, "RME_Q shape");

          // Off-diagonal L branch for Q (m = -2)
          u3::physicalvector p22_0(u3::irrep(2, 2), 0);
          std::vector<int> bk_Q02, kk_Q02;
          auto rme_Q02 = u3::RME_Q(p22_0, p22_2, bk_Q02, kk_Q02);
          require(bk_Q02.size() == 1 && kk_Q02.size() == 2, "RME_Q m=-2 size");
          require(bk_Q02[0] == 0, "RME_Q m=-2 bra_kappa zero-based");
          require(kk_Q02[0] == 0 && kk_Q02[1] == 1, "RME_Q m=-2 ket_kappa zero-based");
        });
  }
} // namespace test
