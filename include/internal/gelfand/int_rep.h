#ifndef U3LIB_INTERNAL_GELFAND_INT_REP_H
#define U3LIB_INTERNAL_GELFAND_INT_REP_H

#include <tuple>

namespace u3
{
  namespace internal
  {
    class gelfand_int_rep
    {
      int n13_, n23_, n33_, n12_, n22_, n11_;

      // helper function for comparison operators:
      auto tie() const { return std::tie(n13_, n23_, n33_, n12_, n22_, n11_); }

    public:
      gelfand_int_rep() : n13_(0), n23_(0), n33_(0), n12_(0), n22_(0), n11_(0)
      { }

      gelfand_int_rep(int n13, int n23, int n33, int n12, int n22, int n11)
          : n13_(n13), n23_(n23), n33_(n33), n12_(n12), n22_(n22), n11_(n11)
      { }

      int n13() const { return n13_; }
      int n23() const { return n23_; }
      int n33() const { return n33_; }
      int n12() const { return n12_; }
      int n22() const { return n22_; }
      int n11() const { return n11_; }

      void n13(int n13) { n13_ = n13; }
      void n23(int n23) { n23_ = n23; }
      void n33(int n33) { n33_ = n33; }
      void n12(int n12) { n12_ = n12; }
      void n22(int n22) { n22_ = n22; }
      void n11(int n11) { n11_ = n11; }

      bool operator==(const gelfand_int_rep& other) const
      {
        return tie() == other.tie();
      }

      bool operator!=(const gelfand_int_rep& other) const
      {
        return !(*this == other);
      }
    };

  } // namespace internal
} // namespace u3

#endif
