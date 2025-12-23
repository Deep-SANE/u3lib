#ifndef U3LIB_IRREP_H
#define U3LIB_IRREP_H

#include <cassert>
#include <cstdint>
#include <istream>
#include <ostream>

#include <qcx_convert.h>
#include <qcx_readwrite.h>

#ifdef HAVE_FMT
#include <fmt/format.h>
#endif

#include "U3gelfand.h"

namespace su3
{
   // Class that represents SU(3) irrep.
  class irrep
  {
  protected:
    // type used to represent labels:
    using rep_type = u3::uqn_t; 

    // SU(3) labels:
    rep_type lambda_, mu_;

  public:
    // default constructor; creates (0,0):
    irrep() : lambda_(0), mu_(0) { }

    irrep(qcx_utils::UInt8 l, qcx_utils::UInt8 m)
        : lambda_(l), mu_(m) // conversion check is implicit
    { }

    // comparison:
    auto operator<=>(const irrep&) const = default;

    // SU(3) getters:
    int lambda() const { return static_cast<int>(lambda_); }
    int lam() const { return lambda(); }
    int mu() const { return static_cast<int>(mu_); }

    // U(3) getters:
    int n13() const { return lambda() + mu(); }
    int n23() const { return mu(); }
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

    // second-order Casimir
    double C2() const
    {
      return (lambda() * lambda() + lambda() * mu() + 
        mu() * mu() + 3.0 * lambda() + 3.0 * mu()) * 2.0 / 3.0;
    }

    // third-order Casimir
    double C3() const
    {
      return (lambda() - mu()) * (lambda() + 2.0 * mu() + 3.0) *
        (2.0 * lambda() + mu() + 3.0) / 9.0;
    }

    // conjugate irrep
    su3::irrep conjugate() const { return su3::irrep(mu(), lambda()); }
    double conjugate_phase() const 
    { 
      return ((lambda() + mu()) % 2 == 0) ? 1.0 : -1.0;
    }

    // write to binary file:
    void writeBinary(std::ostream& s) const
    {
      qcx_utils::writeBinary(s, lambda_);
      qcx_utils::writeBinary(s, mu_);
    }

    // read from binary file:
    void readBinary(std::istream& s)
    {
      qcx_utils::readBinary(s, lambda_);
      qcx_utils::readBinary(s, mu_);
    }
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

    irrep(qcx_utils::UInt8 n1, qcx_utils::UInt8 n2, qcx_utils::UInt8 n3)
        : su3::irrep(n1 - n2, n2 - n3), n33_(n3)
    { 
      assert(n1 >= n2);
      assert(n2 >= n3);
    }

    irrep(qcx_utils::UInt8 l, qcx_utils::UInt8 m) : su3::irrep(l, m), n33_(0) { }

    irrep(const su3::irrep& su3) : su3::irrep(su3), n33_(0) { }

    // SU(3)->U(3) conversion with the total number of quanta provided:
    irrep(const su3::irrep& su3, int N) : irrep(su3)
    {
      assert(N >= 0);
      int n33 = N - su3.lambda() - 2 * su3.mu();
      assert(n33 % 3 == 0);
      n33 /= 3;
      n33_ = n33_ + static_cast<rep_type>(n33);
      assert(quanta() == N);
    }

    // comparison:
    auto operator<=>(const irrep&) const = default;

    // getters:
    int n33() const { return static_cast<int>(n33_); }
    int n13() const { return su3::irrep::n13() + n33(); }
    int n23() const { return su3::irrep::n23() + n33(); }

    // total quanta:
    int quanta() const { return lambda() + 2 * mu() + 3 * n33(); }
    int C1() const { return quanta(); }

    /*
    // check order of quantum numbers
    bool check() { return (n13() >= n23() && n23() >= n33()); }
    */

    // normalize the irrep:
    void normalize() { n33_ = 0; }

    // check for [0,0,0]:
    bool zero() const { return su3::irrep::zero() && (n33_ == 0); }

    // conjugate irrep
    u3::irrep conjugate() const { return u3::irrep(mu(), lambda()); }

    // write to binary file:
    void writeBinary(std::ostream& s) const
    {
      su3::irrep::writeBinary(s);
      qcx_utils::writeBinary(s, n33_);
    }

    // read from binary file:
    void readBinary(std::istream& s)
    {
      su3::irrep::readBinary(s);
      qcx_utils::readBinary(s, n33_);
    }
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
    return fmt::format_to(ctx.out(), "({},{})", i.lambda(), i.mu());
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
    return fmt::format_to(ctx.out(), "[{},{},{}]", i.n13(), i.n23(), i.n33());
  }
};

#endif

#endif
