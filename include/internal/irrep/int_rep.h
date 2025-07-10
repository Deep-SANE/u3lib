#ifndef U3LIB_INTERNAL_IRREP_INT_REP_H
#define U3LIB_INTERNAL_IRREP_INT_REP_H

#include <tuple>

namespace u3
{
  namespace internal
  {
    class irrep_int_rep
    {
      int n13_, n23_, n33_;

      // helper function for comparison operators:
      auto tie() const { return std::tie(n13_, n23_, n33_); }

    public:
      irrep_int_rep() : n13_(0), n23_(0), n33_(0) { }

      irrep_int_rep(int n13, int n23, int n33) : n13_(n13), n23_(n23), n33_(n33)
      { }

      int n13() const { return n13_; }
      int n23() const { return n23_; }
      int n33() const { return n33_; }

      void n13(int n13) { n13_ = n13; }
      void n23(int n23) { n23_ = n23; }
      void n33(int n33) { n33_ = n33; }

      bool operator==(const irrep_int_rep& other) const
      {
        return tie() == other.tie();
      }

      bool operator!=(const irrep_int_rep& other) const
      {
        return !(*this == other);
      }
    };

  } // namespace internal
} // namespace u3

#endif
