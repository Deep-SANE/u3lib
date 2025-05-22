#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace py::literals;

int rhomax(int n31a, int n32a, int n33a, int n31b, int n32b, int n33b, int n31c,
           int n32c, int n33c);

PYBIND11_MODULE(u3lib, m)
{
  m.doc() = "u3lib module";

  m.def("rhomax", &rhomax,
        "Calculate rhomax for U(3) irreps 'a', 'b', and 'c', specified by "
        "their labels",
        "n13a"_a, "n23a"_a, "n33a"_a, "n13b"_a, "n23b"_a, "n33b"_a, "n13c"_a,
        "n23c"_a, "n33c"_a);
}
