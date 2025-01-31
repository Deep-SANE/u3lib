#ifndef U3LIB_GELFAND_GELFAND_H
#define U3LIB_GELFAND_GELFAND_H

#include "int_rep.h"

namespace u3
{
  namespace internal
  {
    template <typename Rep>
    class gelfand : Rep
    {
    public:
      using Rep::n13, Rep::n23, Rep::n33, Rep::n12, Rep::n22, Rep::n11;

      gelfand() = default;

      gelfand(int n13, int n23, int n33, int n12, int n22, int n11)
          : Rep(n13, n23, n33, n12, n22, n11)
      { }

      friend bool operator==(const gelfand& l, const gelfand& r)
      {
         return static_cast<const Rep&>(l) == static_cast<const Rep&>(r);
      }

      friend bool operator!=(const gelfand& l, const gelfand& r)
      {
         return static_cast<const Rep&>(l) != static_cast<const Rep&>(r);
      }

      // check betweenness condition:
      bool check() const
      {
        return n13() >= n12() && n12() >= n23() && n23() >= n22() &&
               n22() >= n33() && n12() >= n11() && n11() >= n22();
      }

      void normalize()
      {
        n13(n13() - n33());
        n23(n23() - n33());
        n12(n12() - n33());
        n22(n22() - n33());
        n11(n11() - n33());
        n33(0);
      }
    };

    using gelfand_int_rep_t = gelfand<gelfand_int_rep>;

  } // namespace internal
} // namespace u3

#endif
