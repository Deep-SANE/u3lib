#pragma once
#include "test_helpers.h"
namespace test
{
  struct TripleChannel
  {
    u3::irrep intermediate;
    int pair_rho, total_rho;
  };

  // Complete multiplicity-resolved basis for (A x B) x C -> total.
  inline auto triple_channels(const u3::irrep& a, const u3::irrep& b,
                              const u3::irrep& c, const u3::irrep& total)
  {
    std::vector<TripleChannel> channels;
    auto [irreps, multiplicities] = u3::decomposition(a, b);
    for (int i = 0; i < int(irreps.size()); ++i)
      for (int r = 0; r < multiplicities[i]; ++r)
        for (int t = 0; t < u3::Rho(irreps[i], c, total); ++t)
        {
          require(irreps[i].lam() <= 5 && irreps[i].mu() <= 5,
                  "Intermediate irrep exceeds test bound");
          channels.push_back({irreps[i], r, t});
        }
    return channels;
  }

  inline void multiplicity_u6z6(Suite& suite, Reference& ref, u3::irrep a,
                                u3::irrep b, u3::irrep c, u3::irrep total,
                                int expected_dimension)
  {
    auto left = triple_channels(a, b, c, total);
    // U6 uses A x (B x C); multiplicity counts agree with (B x C) x A,
    // but coefficient retrieval below preserves the actual coupling order.
    auto rightU = triple_channels(b, c, a, total);
    auto rightZ = triple_channels(a, c, b, total);
    for (bool isZ : {false, true})
    {
      std::string kind = isZ ? "Z6" : "U6";
      suite.run(
          kind + " " + irname(a) + "x" + irname(b) + "x" + irname(c) + " -> " +
              irname(total) +
              ": complete multiplicity matrix against documented values and "
              "orthogonality",
          [&]
          {
            const auto& right = isZ ? rightZ : rightU;
            require(left.size() == size_t(expected_dimension) &&
                        right.size() == left.size(),
                    "Complete recoupling basis dimension");
            Eigen::MatrixXd matrix(left.size(), right.size());
            bool pair_multiplicity = false, total_multiplicity = false;
            for (int i = 0; i < int(left.size()); ++i)
              for (int j = 0; j < int(right.size()); ++j)
              {
                const auto& l = left[i];
                const auto& r = right[j];
                int rAB = u3::Rho(a, b, l.intermediate),
                    rAB_C = u3::Rho(l.intermediate, c, total);
                int rPair = isZ ? u3::Rho(a, c, r.intermediate)
                                : u3::Rho(b, c, r.intermediate);
                int rTotal = isZ ? u3::Rho(r.intermediate, b, total)
                                 : u3::Rho(a, r.intermediate, total);
                pair_multiplicity |= rAB > 1 || rPair > 1;
                total_multiplicity |= rAB_C > 1 || rTotal > 1;
                double value;
                if (isZ)
                {
                  const auto block = u3::u3z6db().get_U3Z6(
                      b, a, total, c, l.intermediate, rAB, rAB_C,
                      r.intermediate, rPair, rTotal);
                  require(block.size() == size_t(rAB * rAB_C * rPair * rTotal),
                          "Z6 block dimensions");
                  value =
                      block(l.pair_rho, l.total_rho, r.pair_rho, r.total_rho);
                }
                else
                {
                  const auto block = u3::u3u6db().get_U3U6(
                      a, b, total, c, l.intermediate, rAB, rAB_C,
                      r.intermediate, rPair, rTotal);
                  require(block.size() == size_t(rAB * rAB_C * rPair * rTotal),
                          "U6 block dimensions");
                  value =
                      block(l.pair_rho, l.total_rho, r.pair_rho, r.total_rho);
                }
                matrix(i, j) = value;
                ref.check(key(kind + "multi", irname(a), irname(b), irname(c),
                              irname(total), irname(l.intermediate), l.pair_rho,
                              l.total_rho, irname(r.intermediate), r.pair_rho,
                              r.total_rho),
                          value);
              }
            require(total_multiplicity,
                    "Expected nontrivial final coupling multiplicity");
            if (a == b && b == c)
              require(pair_multiplicity,
                      "Expected nontrivial pair multiplicity");
            orthogonal(matrix, kind + " complete column orthogonality");
            orthogonal(matrix.transpose(),
                       kind + " complete row orthogonality");
          });
    }
  }

  struct FourChannel
  {
    u3::irrep first, second;
    int first_rho, second_rho, total_rho;
  };
  inline auto four_channels(const u3::irrep& a, const u3::irrep& b,
                            const u3::irrep& c, const u3::irrep& d,
                            const u3::irrep& total)
  {
    std::vector<FourChannel> channels;
    auto [ab, rab] = u3::decomposition(a, b);
    auto [cd, rcd] = u3::decomposition(c, d);
    for (int i = 0; i < int(ab.size()); ++i)
      for (int j = 0; j < int(cd.size()); ++j)
        for (int p = 0; p < rab[i]; ++p)
          for (int q = 0; q < rcd[j]; ++q)
            for (int r = 0; r < u3::Rho(ab[i], cd[j], total); ++r)
            {
              require(ab[i].lam() <= 5 && ab[i].mu() <= 5 && cd[j].lam() <= 5 &&
                          cd[j].mu() <= 5,
                      "Intermediate irrep exceeds test bound");
              channels.push_back({ab[i], cd[j], p, q, r});
            }
    return channels;
  }

  inline void multiplicity_u9(Suite& suite, Reference& ref, u3::irrep total,
                              int expected_dimension)
  {
    u3::irrep f(1, 1);
    suite.run(
        "U9 four (1,1) irreps -> " + irname(total) +
            ": complete multiplicity matrix against documented values and "
            "orthogonality",
        [&]
        {
          auto channels = four_channels(f, f, f, f, total);
          require(channels.size() == size_t(expected_dimension),
                  "Complete U9 basis dimension");
          Eigen::MatrixXd matrix(channels.size(), channels.size());
          bool pair_multiplicity = false, total_multiplicity = false;
          for (int i = 0; i < int(channels.size()); ++i)
            for (int j = 0; j < int(channels.size()); ++j)
            {
              const auto& l = channels[i];
              const auto& r = channels[j];
              int rab = u3::Rho(f, f, l.first), rcd = u3::Rho(f, f, l.second),
                  rl = u3::Rho(l.first, l.second, total);
              int rac = u3::Rho(f, f, r.first), rbd = u3::Rho(f, f, r.second),
                  rr = u3::Rho(r.first, r.second, total);
              pair_multiplicity |= rab > 1 && rcd > 1 && rac > 1 && rbd > 1;
              total_multiplicity |= rl > 1 && rr > 1;
              const auto block = u3::u3u9db().get_U3U9(
                  f, f, f, f, l.first, l.second, r.first, r.second, total, rac,
                  rbd, rl, rab, rcd, rr);
              require(block.size() == size_t(rab * rcd * rl * rac * rbd * rr),
                      "U9 block dimensions");
              // The structured API's axis order is AC, BD, (AB CD), AB, CD, (AC
              // BD).
              double value = block(r.first_rho, r.second_rho, l.total_rho,
                                   l.first_rho, l.second_rho, r.total_rho);
              matrix(i, j) = value;
              ref.check(key("U9multi", irname(f), irname(total),
                            irname(l.first), irname(l.second), l.first_rho,
                            l.second_rho, l.total_rho, irname(r.first),
                            irname(r.second), r.first_rho, r.second_rho,
                            r.total_rho),
                        value);
            }
          require(pair_multiplicity,
                  "Expected multiplicity on all four pair-coupling axes");
          if (total == f)
            require(total_multiplicity,
                    "Expected multiplicity on both final-coupling axes");
          orthogonal(matrix, "U9 complete column orthogonality");
          orthogonal(matrix.transpose(), "U9 complete row orthogonality");
        });
  }

  inline void recoupling(Suite& suite, Reference& ref)
  {
    // Full 2x2 recoupling matrices for three fundamentals with total (1,1).
    // Keep these small analytic anchors alongside the multiplicity cases below.
    for (bool scalar : {true, false})
      suite.run(std::string("U6/Z6 ") +
                    (scalar ? "scalar identity"
                            : "three fundamentals, complete channel matrices") +
                    " against documented values",
                [&]
                {
                  u3::irrep f = scalar ? u3::irrep(0, 0) : u3::irrep(1, 0);
                  u3::irrep total = scalar ? u3::irrep(0, 0) : u3::irrep(1, 1);
                  std::vector<u3::irrep> channels =
                      scalar ? std::vector<u3::irrep>{{0, 0}}
                             : std::vector<u3::irrep>{{2, 0}, {0, 1}};
                  Eigen::MatrixXd U(channels.size(), channels.size()), Z = U;
                  for (int i = 0; i < int(channels.size()); ++i)
                    for (int j = 0; j < int(channels.size()); ++j)
                    {
                      auto u = u3::u3u6db().get(f, f, total, f, channels[i], 1,
                                                1, channels[j], 1, 1);
                      auto z = u3::u3z6db().get(f, f, total, f, channels[i], 1,
                                                1, channels[j], 1, 1);
                      require(u.size() == 1 && z.size() == 1,
                              "Recoupling block dimension");
                      U(i, j) = u[0];
                      Z(i, j) = z[0];
                      ref.check(key("U6", irname(f), irname(total),
                                    irname(channels[i]), irname(channels[j])),
                                u[0]);
                      ref.check(key("Z6", irname(f), irname(total),
                                    irname(channels[i]), irname(channels[j])),
                                z[0]);
                      if (scalar)
                      {
                        near(u[0], 1, "Scalar U6");
                        near(z[0], 1, "Scalar Z6");
                      }
                    }
                  orthogonal(U, "U6 channel orthogonality");
                  orthogonal(Z, "Z6 channel orthogonality");
                });
    for (bool scalar : {true, false})
      suite.run(std::string("U9 ") +
                    (scalar ? "scalar identity"
                            : "four fundamentals, total (0,2), complete "
                              "channel matrix") +
                    " against documented values",
                [&]
                {
                  u3::irrep f = scalar ? u3::irrep(0, 0) : u3::irrep(1, 0);
                  u3::irrep total = scalar ? u3::irrep(0, 0) : u3::irrep(0, 2);
                  std::vector<u3::irrep> channels =
                      scalar ? std::vector<u3::irrep>{{0, 0}}
                             : std::vector<u3::irrep>{{2, 0}, {0, 1}};
                  Eigen::MatrixXd U(channels.size(), channels.size());
                  for (int i = 0; i < int(channels.size()); ++i)
                    for (int j = 0; j < int(channels.size()); ++j)
                    {
                      auto v = u3::u3u9db().get(
                          f, f, f, f, channels[i], channels[i], channels[j],
                          channels[j], total, 1, 1, 1, 1, 1, 1);
                      require(v.size() == 1, "U9 block dimension");
                      U(i, j) = v[0];
                      ref.check(key("U9", irname(f), irname(total),
                                    irname(channels[i]), irname(channels[j])),
                                v[0]);
                      if (scalar)
                        near(v[0], 1, "Scalar U9");
                    }
                  orthogonal(U, "U9 channel orthogonality");
                });
    // Tutorial section 10: these full matrices include its AB=(1,1),
    // BC=(1,2) U6 block and AC=(2,2) Z6 block, with every rho entry.
    multiplicity_u6z6(suite, ref, {2, 0}, {0, 2}, {2, 1}, {2, 1}, 5);
    multiplicity_u6z6(suite, ref, {1, 1}, {1, 1}, {1, 1}, {1, 1}, 8);
    multiplicity_u6z6(suite, ref, {2, 2}, {1, 1}, {1, 1}, {2, 2}, 10);
    multiplicity_u9(suite, ref, {0, 0}, 8);
    multiplicity_u9(suite, ref, {1, 1}, 32);
    suite.run("U9 tutorial scalar-leg example against documented value and "
              "analytic identity",
              [&]
              {
                u3::irrep a(2, 0), zero(0, 0), total(4, 0);
                const auto v = u3::u3u9db().get_U3U9(
                    a, zero, zero, a, a, a, a, a, total, 1, 1, 1, 1, 1, 1);
                require(v.size() == 1, "Tutorial U9 block dimension");
                ref.check("U9tutorial:(2,0):(0,0):(0,0):(2,0):(4,0):",
                          v(0, 0, 0, 0, 0, 0));
                near(v(0, 0, 0, 0, 0, 0), 1, "Tutorial U9 scalar-leg identity");
              });
  }
} // namespace test
