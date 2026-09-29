#include <cassert>

#include <qcx_utils.h>

#include <U3basis.h>
#include <U3decomposition.h>
#include <U3irrep.h>
#include <u3lib.h>
#include <U3U2xU3U2.h>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  u3::irrep u3a(4, 2, 0);
  u3::irrep u3b(4, 2, 0);
  u3::irrep u3c(4, 2, 0);

  auto rhomax = u3::Rho(u3a, u3b, u3c);

  auto gArrA = u3::gelfandBasis(u3a);
  auto gArrB = u3::gelfandBasis(u3a);
  auto gArrC = u3::gelfandBasis(u3a);

  auto coeffs =
      u3::WigU3U2(u3a, u3b, u3c, rhomax, -1, gArrA, gArrB, gArrC);

  assert(coeffs.size() == 111);

  auto acc = coeffs.access(4, 2, 4, 0, 4, 2);
  assert(acc.allowed);
  std::cout << "Single thread:\n" << coeffs.get_row(acc) << std::endl;

  /*
  const auto& [flag, row] =
      u3::wigU3U2db().get(4, 2, 0, 4, 2, 4, 2, 0, 4, 0, 4, 2, 0, 4, 2);
  assert(flag);
  std::cout << row << std::endl;
  */

std::cout << "Multiple threads:" << std::endl;
#pragma omp parallel num_threads(3)
  {
    u3::initialize_thread();

    auto row =
        u3::wigU3U2db().get(4, 2, 0, 4, 2, 4, 2, 0, 4, 0, 4, 2, 0, 4, 2);
    assert(row.size() > 0);
    assert(row.isApprox(coeffs.get_row(acc), 1e-10));
    std::cout << row << std::endl;

    u3::finalize_thread();
  }

  u3::finalize();
}
