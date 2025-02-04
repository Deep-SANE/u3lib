#ifndef U3LIB_IRREP_GELFAND_H
#define U3LIB_IRREP_GELFAND_H

#include <ostream>

#include "int_rep.h"

namespace u3
{
  namespace internal
  {
    /*
     * A class for U(3) irrep [n13,n23,n33]
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
      int lam() const { return n13() - n23(); }
      int mu() const { return n23() - n33(); }

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
    };

    // print u3irrep
    template <typename Rep>
    inline std::ostream& operator<<(std::ostream& s, const irrep<Rep>& u3)
    {
      return s << "[" << u3.n13() << "," << u3.n23() << "," << u3.n33() << "]";
    }

    // type aliases:
    using irrep_int_rep_t = irrep<irrep_int_rep>;

  } // namespace internal
} // namespace u3

#endif
