#ifndef U3LIB_INTERNAL_GELFAND_INT_REP_H
#define U3LIB_INTERNAL_GELFAND_INT_REP_H

#include <tuple>

#include <Int.h>

#include "../../Definitions.h"

namespace u3
{
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

      gelfand_int_rep(utils::UInt8 n13, utils::UInt8 n23, 
                      utils::UInt8 n33, utils::UInt8 n12, 
                      utils::UInt8 n22, utils::UInt8 n11)
          : n13_(n13), n23_(n23), n33_(n33), n12_(n12), n22_(n22), n11_(n11)
      { }

      int n13() const { return static_cast<int>(n13_); }
      int n23() const { return static_cast<int>(n23_); }
      int n33() const { return static_cast<int>(n33_); }
      int n12() const { return static_cast<int>(n12_); }
      int n22() const { return static_cast<int>(n22_); }
      int n11() const { return static_cast<int>(n11_); }

      void n13(utils::UInt8 n13) { n13_ = n13; }
      void n23(utils::UInt8 n23) { n23_ = n23; }
      void n33(utils::UInt8 n33) { n33_ = n33; }
      void n12(utils::UInt8 n12) { n12_ = n12; }
      void n22(utils::UInt8 n22) { n22_ = n22; }
      void n11(utils::UInt8 n11) { n11_ = n11; }

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
