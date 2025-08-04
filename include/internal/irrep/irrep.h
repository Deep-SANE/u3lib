#ifndef U3LIB_INTERNAL_IRREP_IRREP_H
#define U3LIB_INTERNAL_IRREP_IRREP_H

#include <ostream>

#include "int_rep.h"
#include "uint16_t_rep.h"

namespace u3
{
  namespace internal
  {
    /*
     * Class for U(3) irrep [n13,n23,n33]
     *
     * The template argument specifies the internal representation.
     */
    template <typename Rep>
    class irrep : Rep
    {
    public:
      using Rep::n13, Rep::n23, Rep::n33;

      irrep() = default;

      irrep(int n13, int n23, int n33) : Rep(n13, n23, n33) { }

      // construction from Elliott's SU(3) (lambda,mu)
      irrep(int lambda, int mu) : Rep(lambda + mu, mu, 0) { }

      // comparison operators:
      bool operator==(const irrep& other) const
      {
        return static_cast<const Rep&>(*this) == static_cast<const Rep&>(other);
      }

      bool operator!=(const irrep& other) const { return !(*this == other); }

      // lambda, mu
      int lambda() const { return n13() - n23(); }
      int lam() const { return lambda(); }

      int mu() const { return n23() - n33(); }

      // total quanta
      int quanta() const { return n13() + n23() + n33(); }

      // dimension
      int dim() const
      {
        return (1 + n13() - n23()) * (1 + n23() - n33()) * (2 + n13() - n33()) /
               2;
      }

      // check order of quantum numbers
      bool check() { return (n13() >= n23() && n23() >= n33()); }

      // normalize the irrep
      void normalize()
      {
        n13(n13() - n33());
        n23(n23() - n33());
        n33(0);
      }

      bool zero() const { return (n13() == 0) && (n23() == 0) && (n33() == 0); }
    };

    // print u3irrep
    template <typename Rep>
    inline std::ostream& operator<<(std::ostream& s, const irrep<Rep>& u3)
    {
      return s << "[" << u3.n13() << "," << u3.n23() << "," << u3.n33() << "]";
    }

    // type aliases:
    using irrep_int_rep_t = irrep<irrep_int_rep>;
    using irrep_uint16_t_rep_t = irrep<irrep_uint16_t_rep>;

  } // namespace internal
} // namespace u3

namespace su3
{
  namespace internal
  {
    /*
     * Class for SU(3) irrep (lambda,mu)
     *
     * The template argument specifies the internal representation.
     */
    template <typename Rep>
    class irrep : Rep
    {
    public:
      using Rep::lambda, Rep::mu;

      irrep() = default;

      irrep(int lambda, int mu) : Rep(lambda, mu) { }

      template <typename T>
      irrep(const u3::internal::irrep<T>& u3) : Rep(u3.lambda(), u3.mu())
      { }

      // comparison operators:
      bool operator==(const irrep& other) const
      {
        return static_cast<const Rep&>(*this) == static_cast<const Rep&>(other);
      }

      bool operator!=(const irrep& other) const { return !(*this == other); }

      // alternative function name for `lambda`:
      int lam() const { return lambda(); }
      void lam(int lm) { lambda(lm); }

      // U(3) labels:
      int n13() const { return lambda() + mu(); }
      int n23() const { return mu(); }
      int n33() const { return 0; }

      // dimension:
      int dim() const
      {
        return (1 + lambda()) * (1 + mu()) * (2 + lambda() + mu()) / 2;
      }

      bool zero() const { return (lambda() == 0) && (mu() == 0); }
    };

    // type aliases:
    using irrep_uint16_t_rep_t = irrep<irrep_uint16_t_rep>;
  } // namespace internal
} // namespace su3

#endif
