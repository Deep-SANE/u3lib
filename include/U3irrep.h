#ifndef U3LIB_IRREP_H
#define U3LIB_IRREP_H

#include <cassert>
#include <cstdint>
#include <ostream>
#include <tuple>

#include <qcx_convert.h>

#ifdef HAVE_FMT
#include <fmt/format.h>
#endif

namespace su3
{
   // Class that represents SU(3) irrep.
  class irrep
  {
  protected:
    // type used to represent labels:
    using rep_type = uint16_t;

    // SU(3) labels:
    rep_type lambda_, mu_;

  public:
    // default constructor; creates (0,0):
    irrep() : lambda_(0), mu_(0) { }

    irrep(int l, int m)
        : lambda_(qcx_utils::itou16(l)), mu_(qcx_utils::itou16(m))
    { }

    // comparison:
    auto operator<=>(const irrep&) const = default;

    // SU(3) getters:
    int lambda() const { return static_cast<int>(lambda_); }
    int lam() const { return lambda(); }
    int mu() const { return static_cast<int>(mu_); }

    // U(3) getters:
    int n13() const { return lambda_ + mu_; }
    int n23() const { return mu_; }
    int n33() const { return 0; }

    /*
    // setters:
    void lambda(int l) { lambda_ = qcx_utils::itou16(l); }
    void mu(int m) { mu_ = qcx_utils::itou16(m); }
    */

    // irrps's dimesion:
    int dim() const
    {
      return (1 + lambda()) * (1 + mu()) * (2 + lambda() + mu()) / 2;
    }

    // check for (0,0):
    bool zero() const { return (lambda() == 0) && (mu() == 0); }
  };

  // output stream operator:
  inline std::ostream& operator<<(std::ostream& s, const irrep& su3)
  {
    return s << "(" << su3.lambda() << "," << su3.mu() << ")";
  }
} // namespace su3

namespace u3
{
   // Class that represents U(3) irrep.
  class irrep : public su3::irrep
  {
    rep_type n33_;

  public:
    // default constructor; creates [0,0,0]:
    irrep() : n33_(0) { }

    irrep(int n1, int n2, int n3)
        : su3::irrep(n1 - n2, n2 - n3), n33_(qcx_utils::itou16(n3))
    { }

    irrep(int l, int m) : su3::irrep(l, m) { }

    irrep(const su3::irrep& su3) : su3::irrep(su3), n33_(0) { }

    // comparison:
    auto operator<=>(const irrep&) const = default;

    // getters:

    int n13() const { return su3::irrep::n13() + n33_; }
    int n23() const { return su3::irrep::n23() + n33_; }
    int n33() const { return n33_; }

    /*
    // setters:

    void set(int n1, int n2, int n3)
    {
      lambda(n1 - n3);
      mu(n2 - n3);
      n33(n3);
    }

    void n13(int n1) { lambda(n1 - mu_ - n33_); }

    void n23(int n2)
    {
      lambda(n13() - n2);
      mu(n2 - n33_);
    }

    void n33(int n3) { n33_ = qcx_utils::itou16(n3); }
    */

    // total quanta:
    int quanta() const { return lambda_ + 2 * mu_ + 3 * n33_; }

    /*
    // check order of quantum numbers
    bool check() { return (n13() >= n23() && n23() >= n33()); }
    */

    // normalize the irrep:
    void normalize() { n33_ = 0; }

    // check for [0,0,0]:
    bool zero() const { return su3::irrep::zero() && (n33_ == 0); }
  };

  // output stream operator:
  inline std::ostream& operator<<(std::ostream& s, const irrep& u3)
  {
    return s << "[" << u3.n13() << "," << u3.n23() << "," << u3.n33() << "]";
  }
} // namespace u3

#ifdef HAVE_FMT

// formatter for su3::irrep:

template <>
class fmt::formatter<su3::irrep>
{
public:
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  template <typename Context>
  constexpr auto format(const su3::irrep& i, Context& ctx) const
  {
    return format_to(ctx.out(), "({},{})", i.lambda(), i.mu());
  }
};

// formatter for u3::irrep:

template <>
class fmt::formatter<u3::irrep>
{
public:
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  template <typename Context>
  constexpr auto format(const u3::irrep& i, Context& ctx) const
  {
    return format_to(ctx.out(), "[{},{},{}]", i.n13(), i.n23(), i.n33());
  }
};

#endif

#endif
