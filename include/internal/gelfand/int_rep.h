#ifndef U3LIB_INTERNAL_GELFAND_INT_REP_H
#define U3LIB_INTERNAL_GELFAND_INT_REP_H

#include <tuple>

namespace u3
{
  // Define alias datatypes
  using qn_t = int16_t;   // quantum numbers: -32,768 to 32,767
  using mult_t = uint8_t; // multiplicities: 0 to 255

  namespace internal
  {
    class gelfand_int_rep
    {
      qn_t n13_, n23_, n33_, n12_, n22_, n11_;

      // helper function for comparison operators:
      auto tie() const { return std::tie(n13_, n23_, n33_, n12_, n22_, n11_); }

    public:
      gelfand_int_rep() : n13_(0), n23_(0), n33_(0), n12_(0), n22_(0), n11_(0)
      { }

      gelfand_int_rep(qn_t n13, qn_t n23, qn_t n33, qn_t n12, qn_t n22, qn_t n11)
          : n13_(n13), n23_(n23), n33_(n33), n12_(n12), n22_(n22), n11_(n11)
      { }

      qn_t n13() const { return n13_; }
      qn_t n23() const { return n23_; }
      qn_t n33() const { return n33_; }
      qn_t n12() const { return n12_; }
      qn_t n22() const { return n22_; }
      qn_t n11() const { return n11_; }

      void n13(qn_t n13) { n13_ = n13; }
      void n23(qn_t n23) { n23_ = n23; }
      void n33(qn_t n33) { n33_ = n33; }
      void n12(qn_t n12) { n12_ = n12; }
      void n22(qn_t n22) { n22_ = n22; }
      void n11(qn_t n11) { n11_ = n11; }

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
