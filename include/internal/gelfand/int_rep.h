#ifndef U3LIB_INTERNAL_GELFAND_INT_REP_H
#define U3LIB_INTERNAL_GELFAND_INT_REP_H

#include <tuple>

#include <qcx_Int.h>

namespace u3
{
  // Define alias datatypes
  using uqn_t = qcx_utils::UInt8;  // unsigned quantum numbers: 0 to 255
  using qn_t = qcx_utils::Int8;    // signed quantum numbers: -128 to 127
  using mult_t = qcx_utils::UInt8; // multiplicities: 0 to 255

  namespace internal
  {
    class gelfand_int_rep
    {
      uqn_t n13_, n23_, n33_, n12_, n22_, n11_;

      // helper function for comparison operators:
      auto tie() const { return std::tie(n13_, n23_, n33_, n12_, n22_, n11_); }

    public:
      gelfand_int_rep() : n13_(0), n23_(0), n33_(0), n12_(0), n22_(0), n11_(0)
      { }

      gelfand_int_rep(uqn_t n13, uqn_t n23, uqn_t n33, uqn_t n12, uqn_t n22, uqn_t n11)
          : n13_(n13), n23_(n23), n33_(n33), n12_(n12), n22_(n22), n11_(n11)
      { }

      int n13() const { return static_cast<int>(n13_); }
      int n23() const { return static_cast<int>(n23_); }
      int n33() const { return static_cast<int>(n33_); }
      int n12() const { return static_cast<int>(n12_); }
      int n22() const { return static_cast<int>(n22_); }
      int n11() const { return static_cast<int>(n11_); }

      void n13(uqn_t n13) { n13_ = n13; }
      void n23(uqn_t n23) { n23_ = n23; }
      void n33(uqn_t n33) { n33_ = n33; }
      void n12(uqn_t n12) { n12_ = n12; }
      void n22(uqn_t n22) { n22_ = n22; }
      void n11(uqn_t n11) { n11_ = n11; }

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
