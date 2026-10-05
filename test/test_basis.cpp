#include "test_helpers.h"
int main(int argc, char** argv)
{
  u3::initialize(&argc, &argv);
  test::Suite suite;
  for (int la = 0; la <= 5; ++la)
    for (int mu = 0; mu <= 5; ++mu)
    {
      u3::irrep ir(la, mu);
      suite.run(
          "Basis " + test::irname(ir) +
              ": analytic dimension, pattern validity, ordering and SO(3) "
              "branching",
          [&]
          {
            auto gs = u3::gelfandBasis(ir);
            test::require(gs.size() ==
                              size_t((la + 1) * (mu + 1) * (la + mu + 2) / 2),
                          "SU(3) dimension formula");
            int dim = 0;
            for (auto p : u3::physicalBasis(ir))
              dim += (2 * p.L() + 1) * p.kappa();
            test::require(dim == int(gs.size()), "SO(3) branching dimension");
            std::tuple<int, int, int> previous{999, 999, 999};
            for (auto g : gs)
            {
              test::require(g.n13() >= g.n12() && g.n12() >= g.n23() &&
                                g.n23() >= g.n22() && g.n22() >= g.n33() &&
                                g.n12() >= g.n11() && g.n11() >= g.n22(),
                            "Gelfand interlacing");
              std::tuple<int, int, int> current{g.n12(), g.n22(), g.n11()};
              test::require(current < previous,
                            "Deterministic descending Gelfand ordering");
              previous = current;
              u3::gelfand shifted(g.n13() + 1, g.n23() + 1, g.n33() + 1,
                                  g.n12() + 1, g.n22() + 1, g.n11() + 1);
              shifted.normalize();
              test::require(shifted == g, "Determinant-shift normalization");
            }
          });
    }
  u3::finalize();
  return suite.result();
}
