#ifndef U3LIB_PYTHON_WIGU3U2_COEFFS_H
#define U3LIB_PYTHON_WIGU3U2_COEFFS_H

#include <cassert>
#include <optional>

#include <so3Coeffs/SO3_CGC.h>

#include <U3basis.h>
#include <U3decomposition.h>
#include <U3U2xU3U2.h>

#ifdef USE_SLEPC
#include <slepcsys.h>
#endif

class WigU3U2_coeffs
{
  std::optional<u3::WigU3U2_coeffs> coeffs_;

public:
  WigU3U2_coeffs(int n13a, int n23a, int n33a, int n13b, int n23b, int n33b,
                 int n13c, int n23c, int n33c)
  {
#ifdef USE_SLEPC
    SlepcInitializeNoArguments();
#endif

    u3::irrep u3a(n13a, n23a, n33a);
    u3::irrep u3b(n13b, n23b, n33b);
    u3::irrep u3c(n13c, n23c, n33c);

    auto gArrA = u3::gelfandBasis(u3a);
    auto gArrB = u3::gelfandBasis(u3b);
    auto gArrC = u3::gelfandBasis(u3c);

    auto rhomax = u3::Rho(u3a, u3b, u3c);

    assert(rhomax > 0);

    coeffs_.emplace(
        u3::WigU3U2(u3a, u3b, u3c, rhomax, -1, gArrA, gArrB, gArrC));
  }

  auto size() const { return coeffs_->size(); }

  ~WigU3U2_coeffs()
  {
#ifdef USE_SLEPC
    SlepcFinalize();
#endif
  }
};

#endif
