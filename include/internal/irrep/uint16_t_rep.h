#ifndef U3LIB_INTERNAL_IRREP_UINT16_T_REP_H
#define U3LIB_INTERNAL_IRREP_UINT16_T_REP_H

#include <cstdint>
#include <tuple>

#include <qcx_convert.h>

namespace u3
{
  namespace internal
  {
    class irrep_uint16_t_rep
    {
      uint16_t n13_, n23_, n33_;

      // helper function for comparison operators:
      auto tie() const { return std::tie(n13_, n23_, n33_); }

    public:
      irrep_uint16_t_rep() : n13_(0), n23_(0), n33_(0) { }

      irrep_uint16_t_rep(int n13, int n23, int n33)
          : n13_(qcx_utils::itou16(n13)),
            n23_(qcx_utils::itou16(n23)),
            n33_(qcx_utils::itou16(n33))
      { }

      int n13() const { return static_cast<int>(n13_); }
      int n23() const { return static_cast<int>(n23_); }
      int n33() const { return static_cast<int>(n33_); }

      void n13(int n13) { n13_ = qcx_utils::itou16(n13); }
      void n23(int n23) { n23_ = qcx_utils::itou16(n23); }
      void n33(int n33) { n33_ = qcx_utils::itou16(n33); }

      bool operator==(const irrep_uint16_t_rep& other) const
      {
        return tie() == other.tie();
      }

      bool operator!=(const irrep_uint16_t_rep& other) const
      {
        return !(*this == other);
      }
    };

  } // namespace internal
} // namespace u3

namespace su3
{
  namespace internal
  {
    class irrep_uint16_t_rep
    {
      uint16_t lambda_, mu_;

      // helper function for comparison operators:
      auto tie() const { return std::tie(lambda_, mu_); }

    public:
      irrep_uint16_t_rep() : lambda_(0), mu_(0) { }

      irrep_uint16_t_rep(int lambda, int mu)
          : lambda_(qcx_utils::itou16(lambda)), mu_(qcx_utils::itou16(mu))
      { }

      int lambda() const { return static_cast<int>(lambda_); }
      int mu() const { return static_cast<int>(mu_); }

      bool operator==(const irrep_uint16_t_rep& other) const
      {
        return tie() == other.tie();
      }

      bool operator!=(const irrep_uint16_t_rep& other) const
      {
        return !(*this == other);
      }
    };
  } // namespace internal
} // namespace su3

#endif
