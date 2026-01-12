#include <qcx_utils.h>
#include <so3Coeffs/SO3_CGC.h>

#include <u3lib.h>
#include <U3basis.h>
#include <U3decomposition.h>
#include <U3U2xU3U2.h>

int main(int argc, char* argv[])
{
  qcx_utils::initialize(&argc, &argv);
  u3::initialize();

  u3::irrep A(2,0,0), B(2,2,0), C1(4,2,0), C2(0,0,0);

  assert(u3::Rho(A,B,C1) == 1);
  assert(u3::Rho(A,B,C2) == 1);

  std::vector<u3::gelfand> 
    gAs = u3::gelfandBasis(A), gBs = u3::gelfandBasis(B),
    gC1s = u3::gelfandBasis(C1), gC2s = u3::gelfandBasis(C2);

  const auto& gC2 = gC2s[0];
  auto LamC2 = gC2.Lambdax2(), MLamC2 = gC2.MLambdax2();

  assert(gC2s.size() == 1);

  auto& cgcs = qcx_utils::cgcs();
  auto& wigU3U2db = u3::wigU3U2db();

  for (const auto& gC1 : gC1s)
  {
    double sum = 0.0;

    auto LamC1 = gC1.Lambdax2(), MLamC1 = gC1.MLambdax2();

    for (const auto& gA : gAs)
    {
      auto LamA = gA.Lambdax2(), MLamA = gA.MLambdax2();

      for (const auto& gB : gBs)
      {
        auto LamB = gB.Lambdax2(), MLamB = gB.MLambdax2();

        double cgc1 = cgcs.get({LamA, LamB, LamC1, MLamA, MLamB, MLamC1}),
               cgc2 = cgcs.get({LamA, LamB, LamC2, MLamA, MLamB, MLamC2});

        auto [allow1, U3_1] = wigU3U2db.get(2,0,0,gA.n12(),gA.n22(),
                                            2,2,0,gB.n12(),gB.n22(),
                                            4,2,0,gC1.n12(),gC1.n22());

        auto [allow2, U3_2] = wigU3U2db.get(2,0,0,gA.n12(),gA.n22(),
                                            2,2,0,gB.n12(),gB.n22(),
                                            0,0,0,0,0);

        if (!allow1 || !allow2) continue;

        assert(U3_1.size() == 1);
        assert(U3_2.size() == 1);

        sum += (U3_1(0) * cgc1 * U3_2(0) * cgc2);
      }
    }

    std::cout << "gC1: " << gC1 << " => sum = " << sum << std::endl;
  }

  u3::finalize();
  qcx_utils::finalize();
}
