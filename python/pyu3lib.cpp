// ============================================================================
//  pyu3lib.cpp  --  Python bindings for u3lib (pybind11).  Module name: u3lib.
//
//  Design goal: mirror the C++ API 1:1 so users who know u3lib use the same
//  names in Python. No pythonic renaming. Coefficients come back either as
//  numpy arrays (flat, float64) or as the Definitions.h container objects
//  (U3U2/U3SO3/U3U6/U3Z6/U3U9), which give multiplicity-indexed access via
//  __call__ and a flat .data() numpy view.
//
//  ==========================================================================
//  COMPLETE LIST OF EXPOSED NAMES  (Python side; C++ counterpart in u3::/su3::)
//  ==========================================================================
//
//  Module-level functions
//    initialize()                 set up foundation + U(3) DBs (call once first)
//    finalize()                   tear everything down (call once at the end)
//    initialize_thread()          per-thread setup (advanced / multithreaded)
//    finalize_thread()            per-thread teardown
//    clear()                      clear cached coefficients (keep DBs alive)
//    is_initialized()             True only between initialize() and finalize()
//    Rho(la, ma, lb, mb, lc, mc)  coupling multiplicity, by SU(3) labels
//    Rho(a, b, c)                 coupling multiplicity, by irrep objects
//    Kappa(lam, mu, L)            SO(3) branching multiplicity of L in (lam,mu)
//    Kappa(irrep, L)              SO(3) branching multiplicity, by irrep object
//    physicalBasis(irrep)         list of physicalvector (kappa, L) states
//    gelfandBasis(irrep)          list of gelfand states spanning the irrep
//
//  physicalvector(irrep, L)       an SO(3)-reduced basis state
//    irr(), lam(), lambda_(), mu(), n13(), n23(), n33(),
//    kappa(), kappamax(), L(), repr
//
//  Database singletons (the ONLY way to reach coefficients -- no free funcs)
//    wigU3U2db()  -> WigU3U2_db      u3u6db() -> U3U6_db
//    wigU3SO3db() -> WigU3SO3_db     u3z6db() -> U3Z6_db
//                                    u3u9db() -> U3U9_db
//
//  su3.irrep(lambda, mu)          SU(3) irrep       str -> "(lambda,mu)"
//    lambda_(), lam(), mu(), n13(), n23(), n33(), dim(), zero(),
//    C2(), C3(), conjugate(), conjugate_phase(), str, repr,
//    full ordering (== != < <= > >=) and hashable (usable as a dict key)
//
//  irrep(lambda, mu) | irrep(n1, n2, n3)   U(3) irrep; derives from su3.irrep
//    (inherits every su3.irrep method above, and adds / overrides:)
//    n13(), n23(), n33(), quanta(), C1(), normalize(), zero(),
//    conjugate(),  str -> "[n13,n23,n33]",  repr,
//    ordering + hash that compare n33 (not just the su3 base labels)
//
//  gelfand(n13, n23, n33, n12, n22, n11)   Gelfand-Tsetlin state
//    n13(), n23(), n33(), n12(), n22(), n11(), check(), normalize(),
//    lambda_(), mu(), epsilon(), Lambdax2(), MLambdax2(),
//    conjugate(), conjugate_phase(), ==, !=, str, repr, hashable
//    generator matrix elements:
//      e11(), e22(), e33(), e12(), e21(),
//      e23_12(), e23_22(), e32_12(), e32_22(),
//      e13_12(), e13_22(), e31_12(), e31_22()
//
//  Coefficient containers (returned by the get_U3XX accessors below). Each has
//  size(), empty(), data() -> numpy, __getitem__(i) (flat linear index),
//  __call__(...) (multiplicity-indexed), str; plus these shape getters:
//    U3U2   rhomax();                 __call__(rho)
//    U3SO3  kappaAmax(), kappaBmax(), kappaCmax(), rhomax();
//                                     __call__(kappaA, kappaB, kappaC, rho)
//    U3U6   rhoABmax(), rhoAB_Cmax(), rhoBCmax(), rhoA_BCmax();
//                                     __call__(rhoAB, rhoAB_C, rhoBC, rhoA_BC)
//    U3Z6   rhoABmax(), rhoAB_Cmax(), rhoBCmax(), rhoA_BCmax();
//                                     __call__(rhoAB, rhoAB_C, rhoAC, rhoAC_B)
//    U3U9   rhoACmax(), rhoBDmax(), rhoABCDmax(), rhoABmax(), rhoCDmax(),
//           rhoACBDmax();
//                       __call__(rhoAC, rhoBD, rhoABCD, rhoAB, rhoCD, rhoACBD)
//
//  Database methods
//    WigU3U2_db.get(n13a,n23a,n33a,n12a,n22a, ...b, ...c)  -> numpy row over rho
//    WigU3U2_db.get_U3U2(<same 15 labels>)                 -> U3U2 (prints like C++)
//    WigU3SO3_db.get(lam1,mu1,l1, lam2,mu2,l2, lam3,mu3,l3) -> flat numpy
//    WigU3SO3_db.get_U3SO3(<same 9 args>)                  -> U3SO3
//    U3U6_db.get_U3U6(u3A,u3B,u3ABC,u3C, u3AB,rhoABmax,rhoAB_Cmax,
//                     u3BC,rhoBCmax,rhoA_BCmax)             -> U3U6
//    U3Z6_db.get_U3Z6(u3B,u3A,u3ABC,u3C, u3AB,rhoABmax,rhoAB_Cmax,
//                     u3AC,rhoACmax,rhoAC_Bmax)             -> U3Z6
//    U3U9_db.get_U3U9(u3A,u3B,u3C,u3D, u3AB,u3CD, u3AC,u3BD, u3ABCD,
//                     rhoACmax,rhoBDmax,rhoABCDmax, rhoABmax,rhoCDmax,rhoACBDmax)
//                                                           -> U3U9
// ============================================================================

#include <sstream>
#include <stdexcept>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>      // Eigen::RowVectorXd  <-> numpy (wigU3U2db.get)
#include <pybind11/numpy.h>      // py::array_t
#include <pybind11/operators.h>  // py::self == py::self
#include <pybind11/stl.h>        // std::vector <-> list (physicalBasis/gelfandBasis)

#include <u3lib.h>            // umbrella: initialize/finalize (cascades to the
                             // qcx-utils foundation) + all coefficient DBs
#include <U3basis.h>         // Kappa, physicalvector, physical/gelfand bases

namespace py = pybind11;
using namespace py::literals;

namespace
{
  // Copy any contiguous double range (std::vector or std::span) into an owning
  // numpy array. Owning is essential: some DB accessors (U3U9) hand back a
  // std::span view into a thread-local scratch buffer that the next DB call
  // overwrites, so we must copy before returning to Python.
  template <typename Range>
  py::array_t<double> to_numpy(const Range& r)
  {
    return py::array_t<double>(static_cast<py::ssize_t>(std::size(r)),
                               std::data(r));
  }

  // __str__ from a type's C++ operator<<.
  template <typename T>
  std::string to_str(const T& x)
  {
    std::ostringstream os;
    os << x;
    return os.str();
  }

  // Guard for database access. Reaching a coefficient DB before u3::initialize()
  // (or after u3::finalize()) dereferences a null SO3 cache / scratch buffer and
  // segfaults the process -- fatal in a Jupyter kernel. Raise a Python exception
  // instead. u3::initialized_flag() is true only between initialize/finalize.
  void require_init()
  {
    if (!u3::initialized_flag())
      throw std::runtime_error(
          "u3lib is not initialized: call u3lib.initialize() before using a "
          "coefficient database, and do not use the databases after "
          "u3lib.finalize().");
  }
} // namespace

PYBIND11_MODULE(u3lib, m)
{
  m.doc() =
      "u3lib: U(3)/SU(3) coupling and recoupling coefficients.\n"
      "Names mirror the C++ API. Call u3lib.initialize() once before use.";

  // ---- lifecycle -----------------------------------------------------------
  //  u3lib must ALWAYS be initialized before use and finalized after.
  //  initialize() sets up the qcx-utils foundation (Wigner backend tables +
  //  the SO3 CGC/Racah/Fano caches) AND the U(3)-level coefficient databases,
  //  on the calling (main) thread. finalize() tears both down in reverse.
  m.def("initialize", []() { u3::initialize(); },
        "Initialize u3lib (call once before any coefficient call).");
  m.def("finalize", []() { u3::finalize(); },
        "Finalize u3lib (call once when done).");

  // Per-thread setup, for advanced multi-threaded use (mirrors the C++ API).
  m.def("initialize_thread", []() { u3::initialize_thread(); });
  m.def("finalize_thread", []() { u3::finalize_thread(); });

  m.def("clear", &u3::clear, "Clear all cached coefficients (keep DBs alive).");
  m.def("is_initialized", []() { return u3::initialized_flag(); },
        "True between initialize() and finalize(). The coefficient databases "
        "may only be used while this is True.");

  // ==========================================================================
  //  su3::irrep  --  an SU(3) irrep (lambda, mu).  Exposed as u3lib.su3.irrep
  //  so the name mirrors the C++ namespace. u3::irrep derives from it.
  // ==========================================================================
  auto su3 = m.def_submodule("su3", "SU(3) irreps (mirrors the su3:: namespace).");
  py::class_<su3::irrep>(su3, "irrep")
      .def(py::init<>())
      .def(py::init([](int lambda, int mu) { return su3::irrep(lambda, mu); }),
           "lambda"_a, "mu"_a)
      .def("lambda_", &su3::irrep::lambda)   // 'lambda' is a Python keyword
      .def("lam", &su3::irrep::lam)
      .def("mu", &su3::irrep::mu)
      .def("n13", &su3::irrep::n13)
      .def("n23", &su3::irrep::n23)
      .def("n33", &su3::irrep::n33)
      .def("dim", &su3::irrep::dim, "Dimension of the SU(3) irrep.")
      .def("zero", &su3::irrep::zero, "True for the scalar irrep (0,0).")
      .def("C2", &su3::irrep::C2, "Second-order Casimir.")
      .def("C3", &su3::irrep::C3, "Third-order Casimir.")
      .def("conjugate", &su3::irrep::conjugate)
      .def("conjugate_phase", &su3::irrep::conjugate_phase)
      .def(py::self == py::self)   // full ordering, from the C++ operator<=>
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self <= py::self)
      .def(py::self > py::self)
      .def(py::self >= py::self)
      .def("__hash__",
           [](const su3::irrep& i) {
             return py::hash(py::make_tuple(i.lam(), i.mu()));
           })
      .def("__str__", &to_str<su3::irrep>)     // "(lambda,mu)"
      .def("__repr__", [](const su3::irrep& i) {
        return "su3.irrep(" + std::to_string(i.lam()) + "," +
               std::to_string(i.mu()) + ")";
      });

  // ==========================================================================
  //  u3::irrep  --  a U(3) irrep. Labelled either by SU(3) content (lambda,mu)
  //  with n33 = 0, or by the three U(3) row labels (n1 >= n2 >= n3).
  // ==========================================================================
  py::class_<u3::irrep, su3::irrep>(m, "irrep")
      .def(py::init<>())
      .def(py::init([](int lambda, int mu) { return u3::irrep(lambda, mu); }),
           "lambda"_a, "mu"_a, "U(3) irrep from SU(3) labels (n33 = 0).")
      .def(py::init([](int n1, int n2, int n3) { return u3::irrep(n1, n2, n3); }),
           "n1"_a, "n2"_a, "n3"_a, "U(3) irrep from row labels n1 >= n2 >= n3.")
      // n13/n23/n33 shadow the su3 base (they include n33 here):
      .def("n13", &u3::irrep::n13)
      .def("n23", &u3::irrep::n23)
      .def("n33", &u3::irrep::n33)
      .def("quanta", &u3::irrep::quanta, "Total number of oscillator quanta.")
      .def("C1", &u3::irrep::C1, "First-order Casimir (= total quanta).")
      .def("normalize", &u3::irrep::normalize, "Set n33 to 0 in place.")
      .def("zero", &u3::irrep::zero, "True for the scalar irrep [0,0,0].")
      .def("conjugate", &u3::irrep::conjugate)
      // u3-level comparisons (from u3::irrep::operator<=>): these compare n33,
      // unlike the inherited su3 == which only sees (lambda, mu).
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self <= py::self)
      .def(py::self > py::self)
      .def(py::self >= py::self)
      .def("__hash__",
           [](const u3::irrep& i) {
             return py::hash(py::make_tuple(i.n13(), i.n23(), i.n33()));
           })
      .def("__str__", &to_str<u3::irrep>)      // "[n13,n23,n33]"
      .def("__repr__", [](const u3::irrep& i) {
        return "irrep(" + std::to_string(i.n13()) + "," +
               std::to_string(i.n23()) + "," + std::to_string(i.n33()) + ")";
      });

  // ==========================================================================
  //  u3::gelfand  --  a Gelfand basis state
  //        | n13 n23 n33 >
  //        |   n12 n22   >
  //        |     n11     >
  //  Carries the SU(2) subgroup labels and the U(3) generator matrix elements.
  // ==========================================================================
  py::class_<u3::gelfand>(m, "gelfand")
      .def(py::init([](int n13, int n23, int n33, int n12, int n22, int n11) {
             return u3::gelfand(n13, n23, n33, n12, n22, n11);
           }),
           "n13"_a, "n23"_a, "n33"_a, "n12"_a, "n22"_a, "n11"_a)
      // labels (getters; the C++ names are overloaded getter/setter, so bind
      // via lambdas to pick the const getter):
      .def("n13", [](const u3::gelfand& g) { return g.n13(); })
      .def("n23", [](const u3::gelfand& g) { return g.n23(); })
      .def("n33", [](const u3::gelfand& g) { return g.n33(); })
      .def("n12", [](const u3::gelfand& g) { return g.n12(); })
      .def("n22", [](const u3::gelfand& g) { return g.n22(); })
      .def("n11", [](const u3::gelfand& g) { return g.n11(); })
      .def("check", &u3::gelfand::check, "Betweenness (validity) condition.")
      .def("normalize", &u3::gelfand::normalize)
      .def("lambda_", &u3::gelfand::lambda)
      .def("lam", &u3::gelfand::lambda)
      .def("mu", &u3::gelfand::mu)
      .def("epsilon", &u3::gelfand::epsilon, "U(1) weight epsilon.")
      .def("Lambdax2", &u3::gelfand::Lambdax2, "2 * SU(2) angular momentum.")
      .def("MLambdax2", &u3::gelfand::MLambdax2, "2 * SU(2) projection.")
      .def("conjugate", &u3::gelfand::conjugate)
      .def("conjugate_phase", &u3::gelfand::conjugate_phase)
      // U(3) generator matrix elements on this state:
      .def("e11", &u3::gelfand::e11)
      .def("e22", &u3::gelfand::e22)
      .def("e33", &u3::gelfand::e33)
      .def("e12", &u3::gelfand::e12)
      .def("e21", &u3::gelfand::e21)
      .def("e23_12", &u3::gelfand::e23_12)
      .def("e23_22", &u3::gelfand::e23_22)
      .def("e32_12", &u3::gelfand::e32_12)
      .def("e32_22", &u3::gelfand::e32_22)
      .def("e13_12", &u3::gelfand::e13_12)
      .def("e13_22", &u3::gelfand::e13_22)
      .def("e31_12", &u3::gelfand::e31_12)
      .def("e31_22", &u3::gelfand::e31_22)
      .def(py::self == py::self)   // gelfand has == / != only (no C++ <=>)
      .def(py::self != py::self)
      .def("__hash__",
           [](const u3::gelfand& g) {
             return py::hash(py::make_tuple(g.n13(), g.n23(), g.n33(), g.n12(),
                                            g.n22(), g.n11()));
           })
      .def("__str__", &to_str<u3::gelfand>)
      .def("__repr__", [](const u3::gelfand& g) {
        return "gelfand(" + to_str(g) + ")";
      });

  // ---- Rho : coupling multiplicity of (la,ma) x (lb,mb) -> (lc,mc) ----------
  m.def("Rho",
        [](int la, int ma, int lb, int mb, int lc, int mc) {
          return u3::Rho(u3::irrep(la, ma), u3::irrep(lb, mb),
                         u3::irrep(lc, mc));
        },
        "la"_a, "ma"_a, "lb"_a, "mb"_a, "lc"_a, "mc"_a,
        "Coupling multiplicity rho (by SU(3) labels).");
  m.def("Rho",
        [](const u3::irrep& a, const u3::irrep& b, const u3::irrep& c) {
          return u3::Rho(a, b, c);
        },
        "a"_a, "b"_a, "c"_a, "Coupling multiplicity rho (by irrep objects).");

  // ---- Kappa : SO(3) branching multiplicity of L in (lambda, mu) -----------
  //  (kappa_max: how many times angular momentum L appears in the irrep).
  m.def("Kappa",
        [](int lam, int mu, int L) { return u3::Kappa(u3::irrep(lam, mu), L); },
        "lam"_a, "mu"_a, "L"_a,
        "SO(3) branching multiplicity kappa_max of L in (lambda,mu).");
  m.def("Kappa", [](const u3::irrep& a, int L) { return u3::Kappa(a, L); },
        "irrep"_a, "L"_a, "SO(3) branching multiplicity (by irrep object).");

  // ---- physicalvector : an SO(3)-reduced basis state |(lambda,mu) kappa L> --
  py::class_<u3::physicalvector>(m, "physicalvector")
      .def(py::init([](const u3::irrep& irr, int L) {
             return u3::physicalvector(irr, L);
           }),
           "irrep"_a, "L"_a)
      .def("irr", &u3::physicalvector::irr,
           py::return_value_policy::reference_internal, "The U(3) irrep.")
      .def("lam", &u3::physicalvector::lam)
      .def("lambda_", &u3::physicalvector::lambda)
      .def("mu", &u3::physicalvector::mu)
      .def("n13", &u3::physicalvector::n13)
      .def("n23", &u3::physicalvector::n23)
      .def("n33", &u3::physicalvector::n33)
      .def("kappa", &u3::physicalvector::kappa, "SO(3) multiplicity of L.")
      .def("kappamax", &u3::physicalvector::kappamax)
      .def("L", &u3::physicalvector::L)
      .def("__repr__", [](const u3::physicalvector& p) {
        return "physicalvector((" + std::to_string(p.lam()) + "," +
               std::to_string(p.mu()) + "), kappa=" + std::to_string(p.kappa()) +
               ", L=" + std::to_string(p.L()) + ")";
      });

  // ---- basis enumerators (return Python lists) -----------------------------
  m.def("physicalBasis", &u3::physicalBasis, "irrep"_a,
        "List of physicalvector states (kappa, L) spanning the irrep.");
  m.def("gelfandBasis", &u3::gelfandBasis, "irrep"_a,
        "List of gelfand states spanning the irrep.");

  // ==========================================================================
  //  Coefficient container classes (Definitions.h). These are what the
  //  structured get_U3XX() accessors return: they own a flat coefficient
  //  buffer and expose multiplicity-indexed access via __call__.
  //  __getitem__ is flat (linear index); .data() is the whole buffer as numpy.
  // ==========================================================================
  // U3U2: a 1-D coefficient array (a single multiplicity rho). Included for
  // consistency with the family; __call__ and __getitem__ both index by rho.
  using U3U2_t = u3::U3U2<u3::coeffs_t>;
  py::class_<U3U2_t>(m, "U3U2")
      .def("size", &U3U2_t::size)
      .def("empty", &U3U2_t::empty)
      .def("rhomax", &U3U2_t::rhomax)
      .def("__call__", [](const U3U2_t& o, size_t rho) { return o(rho); },
           "rho"_a)
      .def("__getitem__", [](const U3U2_t& o, size_t i) { return o[i]; })
      .def("data", [](const U3U2_t& o) { return to_numpy(o.data()); })
      .def("__str__", &to_str<U3U2_t>);

  using U3SO3_t = u3::U3SO3<u3::coeffs_t>;
  py::class_<U3SO3_t>(m, "U3SO3")
      .def("size", &U3SO3_t::size)
      .def("empty", &U3SO3_t::empty)
      .def("kappaAmax", &U3SO3_t::kappaAmax)
      .def("kappaBmax", &U3SO3_t::kappaBmax)
      .def("kappaCmax", &U3SO3_t::kappaCmax)
      .def("rhomax", &U3SO3_t::rhomax)
      .def("__call__",
           [](const U3SO3_t& o, size_t k1, size_t k2, size_t k3, size_t rho) {
             return o(k1, k2, k3, rho);
           },
           "kappaA"_a, "kappaB"_a, "kappaC"_a, "rho"_a)
      .def("__getitem__", [](const U3SO3_t& o, size_t i) { return o[i]; })
      .def("data", [](const U3SO3_t& o) { return to_numpy(o.data()); })
      .def("__str__", &to_str<U3SO3_t>);

  using U3U6_t = u3::U3U6<u3::coeffs_t>;
  py::class_<U3U6_t>(m, "U3U6")
      .def("size", &U3U6_t::size)
      .def("empty", &U3U6_t::empty)
      .def("rhoABmax", &U3U6_t::rhoABmax)
      .def("rhoAB_Cmax", &U3U6_t::rhoAB_Cmax)
      .def("rhoBCmax", &U3U6_t::rhoBCmax)
      .def("rhoA_BCmax", &U3U6_t::rhoA_BCmax)
      .def("__call__",
           [](const U3U6_t& o, size_t rhoAB, size_t rhoAB_C, size_t rhoBC,
              size_t rhoA_BC) { return o(rhoAB, rhoAB_C, rhoBC, rhoA_BC); },
           "rhoAB"_a, "rhoAB_C"_a, "rhoBC"_a, "rhoA_BC"_a)
      .def("__getitem__", [](const U3U6_t& o, size_t i) { return o[i]; })
      .def("data", [](const U3U6_t& o) { return to_numpy(o.data()); })
      .def("__str__", &to_str<U3U6_t>);

  using U3Z6_t = u3::U3Z6<u3::coeffs_t>;
  py::class_<U3Z6_t>(m, "U3Z6")
      .def("size", &U3Z6_t::size)
      .def("empty", &U3Z6_t::empty)
      .def("rhoABmax", &U3Z6_t::rhoABmax)
      .def("rhoAB_Cmax", &U3Z6_t::rhoAB_Cmax)
      .def("rhoBCmax", &U3Z6_t::rhoBCmax)     // header getter name (dim 2 = rhoAC)
      .def("rhoA_BCmax", &U3Z6_t::rhoA_BCmax)
      .def("__call__",
           [](const U3Z6_t& o, size_t rhoAB, size_t rhoAB_C, size_t rhoAC,
              size_t rhoAC_B) { return o(rhoAB, rhoAB_C, rhoAC, rhoAC_B); },
           "rhoAB"_a, "rhoAB_C"_a, "rhoAC"_a, "rhoAC_B"_a)
      .def("__getitem__", [](const U3Z6_t& o, size_t i) { return o[i]; })
      .def("data", [](const U3Z6_t& o) { return to_numpy(o.data()); })
      .def("__str__", &to_str<U3Z6_t>);

  using U3U9_t = u3::U3U9<u3::coeffs_t>;
  py::class_<U3U9_t>(m, "U3U9")
      .def("size", &U3U9_t::size)
      .def("empty", &U3U9_t::empty)
      .def("rhoACmax", &U3U9_t::rhoACmax)
      .def("rhoBDmax", &U3U9_t::rhoBDmax)
      .def("rhoABCDmax", &U3U9_t::rhoABCDmax)
      .def("rhoABmax", &U3U9_t::rhoABmax)
      .def("rhoCDmax", &U3U9_t::rhoCDmax)
      .def("rhoACBDmax", &U3U9_t::rhoACBDmax)
      .def("__call__",
           [](const U3U9_t& o, size_t rhoAC, size_t rhoBD, size_t rhoABCD,
              size_t rhoAB, size_t rhoCD, size_t rhoACBD) {
             return o(rhoAC, rhoBD, rhoABCD, rhoAB, rhoCD, rhoACBD);
           },
           "rhoAC"_a, "rhoBD"_a, "rhoABCD"_a, "rhoAB"_a, "rhoCD"_a, "rhoACBD"_a)
      .def("__getitem__", [](const U3U9_t& o, size_t i) { return o[i]; })
      .def("data", [](const U3U9_t& o) { return to_numpy(o.data()); })
      .def("__str__", &to_str<U3U9_t>);

  // ==========================================================================
  //  Coefficient databases (singletons). All coefficients are reached through
  //  these -- there are no standalone free functions in the public API.
  // ==========================================================================

  // ---- U(3) > U(2) reduced coupling coefficients ---------------------------
  //  get(...) takes, per irrep a/b/c: (n13, n23, n33, n12, n22).
  //  Returns the coefficient row over rho as a numpy array (1-D).
  py::class_<u3::WigU3U2_db, std::unique_ptr<u3::WigU3U2_db, py::nodelete>>(
      m, "WigU3U2_db")
      .def("get",
           [](u3::WigU3U2_db& db, int n13a, int n23a, int n33a, int n12a,
              int n22a, int n13b, int n23b, int n33b, int n12b, int n22b,
              int n13c, int n23c, int n33c, int n12c, int n22c) {
             require_init();
             return db.get(n13a, n23a, n33a, n12a, n22a, n13b, n23b, n33b, n12b,
                           n22b, n13c, n23c, n33c, n12c, n22c);
           },
           "n13a"_a, "n23a"_a, "n33a"_a, "n12a"_a, "n22a"_a,
           "n13b"_a, "n23b"_a, "n33b"_a, "n12b"_a, "n22b"_a,
           "n13c"_a, "n23c"_a, "n33c"_a, "n12c"_a, "n22c"_a,
           "Coefficient row over rho (numpy array).")
      // Structured U3U2 object: same coefficients over rho, with .data() and a
      // __str__/print that reproduce the C++ U3U2 layout.
      .def("get_U3U2",
           [](u3::WigU3U2_db& db, int n13a, int n23a, int n33a, int n12a,
              int n22a, int n13b, int n23b, int n33b, int n12b, int n22b,
              int n13c, int n23c, int n33c, int n12c, int n22c) {
             require_init();
             Eigen::RowVectorXd row =
                 db.get(n13a, n23a, n33a, n12a, n22a, n13b, n23b, n33b, n12b,
                        n22b, n13c, n23c, n33c, n12c, n22c);
             u3::coeffs_t v(row.data(), row.data() + row.size());
             return u3::U3U2<u3::coeffs_t>(std::move(v));
           },
           "n13a"_a, "n23a"_a, "n33a"_a, "n12a"_a, "n22a"_a,
           "n13b"_a, "n23b"_a, "n33b"_a, "n12b"_a, "n22b"_a,
           "n13c"_a, "n23c"_a, "n33c"_a, "n12c"_a, "n22c"_a,
           "Structured U3U2 object (over rho; prints like C++).");
  m.def("wigU3U2db", &u3::wigU3U2db, py::return_value_policy::reference,
        "The global U(3)>U(2) coefficient database (singleton).");

  // ---- U(3) > SO(3) reduced coupling coefficients --------------------------
  //  get(...) per irrep: (lambda, mu, L). Flat numpy row, or the structured
  //  U3SO3 object via get_U3SO3(...) (indexed by kappaA,kappaB,kappaC,rho).
  py::class_<u3::WigU3SO3_db, std::unique_ptr<u3::WigU3SO3_db, py::nodelete>>(
      m, "WigU3SO3_db")
      .def(
          "get",
          [](u3::WigU3SO3_db& db, int lam1, int mu1, int l1, int lam2, int mu2,
             int l2, int lam3, int mu3, int l3) {
            require_init();
            return to_numpy(db.get(lam1, mu1, l1, lam2, mu2, l2, lam3, mu3, l3));
          },
          "lam1"_a, "mu1"_a, "l1"_a, "lam2"_a, "mu2"_a, "l2"_a, "lam3"_a,
          "mu3"_a, "l3"_a, "Coefficients as a flat numpy array.")
      .def("get_U3SO3",
           [](u3::WigU3SO3_db& db, int lam1, int mu1, int l1, int lam2, int mu2,
              int l2, int lam3, int mu3, int l3) {
             require_init();
             return db.get_U3SO3(lam1, mu1, l1, lam2, mu2, l2, lam3, mu3, l3);
           },
           "lam1"_a, "mu1"_a, "l1"_a, "lam2"_a, "mu2"_a, "l2"_a, "lam3"_a,
           "mu3"_a, "l3"_a, "Structured U3SO3 object (kappaA,kappaB,kappaC,rho).");
  m.def("wigU3SO3db", &u3::wigU3SO3db, py::return_value_policy::reference,
        "The global U(3)>SO(3) coefficient database (singleton).");

  // ---- U(3) 6-(lambda,mu) recoupling: U-coefficients (U3U6) ----------------
  //  get_U3U6(...) needs the four coupling multiplicities (from Rho), matching
  //  the C++ signature. Returns a U3U6 object.
  py::class_<u3::U3U6_db, std::unique_ptr<u3::U3U6_db, py::nodelete>>(
      m, "U3U6_db")
      .def("get_U3U6",
           [](u3::U3U6_db& db, const u3::irrep& u3A, const u3::irrep& u3B,
              const u3::irrep& u3ABC, const u3::irrep& u3C, const u3::irrep& u3AB,
              int rhoABmax, int rhoAB_Cmax, const u3::irrep& u3BC, int rhoBCmax,
              int rhoA_BCmax) {
             require_init();
             return db.get_U3U6(u3A, u3B, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                u3BC, rhoBCmax, rhoA_BCmax);
           },
           "u3A"_a, "u3B"_a, "u3ABC"_a, "u3C"_a,
           "u3AB"_a, "rhoABmax"_a, "rhoAB_Cmax"_a,
           "u3BC"_a, "rhoBCmax"_a, "rhoA_BCmax"_a,
           "Structured U3U6 recoupling object.");
  m.def("u3u6db", &u3::u3u6db, py::return_value_policy::reference,
        "The global U(3) U-coefficient database (singleton).");

  // ---- U(3) 6-(lambda,mu) recoupling: Z-coefficients (U3Z6) ----------------
  py::class_<u3::U3Z6_db, std::unique_ptr<u3::U3Z6_db, py::nodelete>>(
      m, "U3Z6_db")
      .def("get_U3Z6",
           [](u3::U3Z6_db& db, const u3::irrep& u3B, const u3::irrep& u3A,
              const u3::irrep& u3ABC, const u3::irrep& u3C, const u3::irrep& u3AB,
              int rhoABmax, int rhoAB_Cmax, const u3::irrep& u3AC, int rhoACmax,
              int rhoAC_Bmax) {
             require_init();
             return db.get_U3Z6(u3B, u3A, u3ABC, u3C, u3AB, rhoABmax, rhoAB_Cmax,
                                u3AC, rhoACmax, rhoAC_Bmax);
           },
           "u3B"_a, "u3A"_a, "u3ABC"_a, "u3C"_a,
           "u3AB"_a, "rhoABmax"_a, "rhoAB_Cmax"_a,
           "u3AC"_a, "rhoACmax"_a, "rhoAC_Bmax"_a,
           "Structured U3Z6 recoupling object.");
  m.def("u3z6db", &u3::u3z6db, py::return_value_policy::reference,
        "The global U(3) Z-coefficient database (singleton).");

  // ---- U(3) 9-(lambda,mu) recoupling coefficients (U3U9) -------------------
  //  The C++ get_U3U9 returns a view (std::span) into a scratch buffer; we copy
  //  it into an owning U3U9 so the Python object stays valid across DB calls.
  py::class_<u3::U3U9_db, std::unique_ptr<u3::U3U9_db, py::nodelete>>(
      m, "U3U9_db")
      .def(
          "get_U3U9",
          [](u3::U3U9_db& db, const u3::irrep& u3A, const u3::irrep& u3B,
             const u3::irrep& u3C, const u3::irrep& u3D, const u3::irrep& u3AB,
             const u3::irrep& u3CD, const u3::irrep& u3AC, const u3::irrep& u3BD,
             const u3::irrep& u3ABCD, int rhoACmax, int rhoBDmax, int rhoABCDmax,
             int rhoABmax, int rhoCDmax, int rhoACBDmax) {
            require_init();
            auto src =
                db.get_U3U9(u3A, u3B, u3C, u3D, u3AB, u3CD, u3AC, u3BD, u3ABCD,
                            rhoACmax, rhoBDmax, rhoABCDmax, rhoABmax, rhoCDmax,
                            rhoACBDmax);
            // Copy the (possibly span-backed) buffer into an owning container.
            u3::coeffs_t buf(src.data().begin(), src.data().end());
            return U3U9_t(std::move(buf), src.rhoACmax(), src.rhoBDmax(),
                          src.rhoABCDmax(), src.rhoABmax(), src.rhoCDmax(),
                          src.rhoACBDmax());
          },
          "u3A"_a, "u3B"_a, "u3C"_a, "u3D"_a, "u3AB"_a, "u3CD"_a, "u3AC"_a,
          "u3BD"_a, "u3ABCD"_a, "rhoACmax"_a, "rhoBDmax"_a, "rhoABCDmax"_a,
          "rhoABmax"_a, "rhoCDmax"_a, "rhoACBDmax"_a,
          "Structured U3U9 recoupling object (owning copy).");
  m.def("u3u9db", &u3::u3u9db, py::return_value_policy::reference,
        "The global U(3) 9-(lambda,mu) coefficient database (singleton).");
}
