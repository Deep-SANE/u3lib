#pragma once
#include "test_helpers.h"
namespace test
{
  inline void conjugation(Suite& suite, Reference& ref)
  {
    for (int la = 0; la <= 5; ++la)
      for (int mu = 0; mu <= 5; ++mu)
      {
        u3::irrep ir(la, mu), bar(mu, la);
        suite.run(
            "Conjugation " + irname(ir) +
                ": documented Gelfand signs, physical action, all L/M and "
                "lowering",
            [&]
            {
              auto gs = u3::gelfandBasis(ir), bs = u3::gelfandBasis(bar);
              std::vector<int> permutation(gs.size());
              for (int i = 0; i < int(gs.size()); ++i)
              {
                auto g = gs[i], h = g.conjugate();
                require(h.conjugate() == g, "Gelfand double conjugation");
                near(g.conjugate_phase() * h.conjugate_phase(), 1,
                     "Gelfand phase squared");
                ref.check(key(irname(ir), "GT", int(g.n12()), int(g.n22()),
                              int(g.n11())),
                          g.conjugate_phase());
                auto it = std::find(bs.begin(), bs.end(), h);
                require(it != bs.end(), "Missing conjugate pattern");
                permutation[i] = std::distance(bs.begin(), it);
              }
              for (auto p : u3::physicalBasis(ir))
              {
                auto [ids, vs] = u3::phys_to_cano_allM(p, gs);
                auto [bids, bvs] =
                    u3::phys_to_cano_allM(u3::physicalvector(bar, p.L()), bs);
                auto lowerids = ids[0];
                Eigen::MatrixXd lowerv = vs[0], D;
                for (int M = p.L(); M >= -p.L(); --M)
                {
                  auto f = full(ids[p.L() - M], vs[p.L() - M], gs.size());
                  auto bf = full(bids[p.L() + M], bvs[p.L() + M], bs.size());
                  auto [singleids, singlev] = u3::phys_to_cano(p, M, gs);
                  near((f - full(singleids, singlev, gs.size())).norm(), 0,
                       "Single/all-M agreement", 1e-9, 0);
                  orthogonal(f, "Physical basis orthonormality");
                  if (M < p.L())
                  {
                    auto next = u3::lowerM(p, M, gs, lowerids, lowerv);
                    lowerids = next.first;
                    lowerv = next.second;
                  }
                  near((f - full(lowerids, lowerv, gs.size())).norm(), 0,
                       "Explicit L_- lowering", 1e-9, 0);
                  Eigen::MatrixXd mapped =
                      Eigen::MatrixXd::Zero(bs.size(), f.cols());
                  for (int i = 0; i < int(gs.size()); ++i)
                    mapped.row(permutation[i]) =
                        gs[i].conjugate_phase() * f.row(i);
                  double phase = sign(la + mu - p.L() + M);
                  if (la != mu)
                  {
                    near((mapped - phase * bf).norm(), 0,
                         "Distinct-irrep conjugation", 1e-9, 0);
                    if (M == p.L())
                      ref.check(key(irname(ir), "phase", p.L()), phase);
                  }
                  else
                  {
                    Eigen::MatrixXd d = phase * bf.transpose() * mapped;
                    if (M == p.L())
                    {
                      D = d;
                      for (int a = 0; a < d.rows(); ++a)
                        for (int b = 0; b < d.cols(); ++b)
                          ref.check(key(irname(ir), "D", p.L(), a, b), d(a, b));
                    }
                    near((d - D).norm(), 0, "Conjugation matrix M independence",
                         1e-9, 0);
                    near((mapped - phase * bf * D).norm(), 0,
                         "Conjugation reconstruction", 1e-9, 0);
                    near((D - D.transpose()).norm(), 0, "Conjugation symmetry",
                         1e-9, 0);
                    near((D * D - Eigen::MatrixXd::Identity(D.rows(), D.cols()))
                             .norm(),
                         0, "Conjugation squared", 1e-9, 0);
                  }
                }
                if (la == 4 && mu == 4 && p.L() == 2)
                {
                  near(D(0, 0), 4.0 / 7, "(4,4) L=2 diagonal");
                  near(D(0, 1), std::sqrt(33.0) / 7, "(4,4) L=2 mixing");
                  near(D(1, 0), std::sqrt(33.0) / 7, "(4,4) L=2 mixing");
                  near(D(1, 1), -4.0 / 7, "(4,4) L=2 diagonal");
                }
              }
            });
      }
  }
} // namespace test
