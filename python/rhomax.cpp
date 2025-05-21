#include <cassert>

#include <U3decomposition.h>

#include <pybind11/pybind11.h>

//extern "C"
//{

int rhomax(int n31a, int n32a, int n33a, int n31b, int n32b, int n33b, int n31c,
           int n32c, int n33c)
{
  u3::irrep u3a(n31a, n32a, n33a);
  u3::irrep u3b(n31b, n32b, n33b);
  u3::irrep u3c(n31c, n32c, n33c);

  assert(u3a.check());
  assert(u3b.check());
  assert(u3c.check());

  return u3::Rho(u3a, u3b, u3c);
}

//}

PYBIND11_MODULE(u3lib, m) { m.def("rhomax", &rhomax); }
