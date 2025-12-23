#ifndef U3LIB_U3BASIS_H
#define U3LIB_U3BASIS_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "U3gelfand.h"
#include "U3irrep.h"

namespace su3
{
  inline int kmax(int lm, int mu, int L)
  {
    assert((lm + mu - L) >= 0);

    return (lm + mu - L) / 2 - std::max(0, (lm + 1 - L) / 2) -
           std::max(0, (mu + 1 - L) / 2) + 1;
  }

  inline int kmax(const su3::irrep& ir, int L)
  {
    return su3::kmax(ir.lambda(), ir.mu(), L);
  }

  // alternatives called 'Kappa':
  inline int Kappa(int lm, int mu, int L) 
  { 
    return su3::kmax(lm, mu, L); 
  }

  inline int Kappa(const su3::irrep& ir, int L) 
  { 
    return su3::kmax(ir, L); 
  }

  // alternatives called 'mult':
  inline int mult(int lm, int mu, int L) 
  { 
    return su3::kmax(lm, mu, L); 
  }

  inline int mult(const su3::irrep& ir, int L) 
  { 
    return su3::kmax(ir, L); 
  }
} // namespace su3

namespace u3
{
  // compare two gelfand patterns of the same irrep
  inline bool E11(const u3::gelfand& i, const u3::gelfand& f)
  {
    return i.n12() == f.n12() && i.n22() == f.n22() && i.n11() == f.n11();
  }

  // define action of the raising generators, that of the lowering generators
  // is simply swapping the initial and final states
  inline bool E12(const u3::gelfand& i, const u3::gelfand& f)
  {
    return i.n12() == f.n12() && i.n22() == f.n22() && (i.n11() + 1) == f.n11();
  }

  inline bool E23_12(const u3::gelfand& i, const u3::gelfand& f)
  {
    return (i.n12() + 1) == f.n12() && i.n22() == f.n22() && i.n11() == f.n11();
  }

  inline bool E23_22(const u3::gelfand& i, const u3::gelfand& f)
  {
    return i.n12() == f.n12() && (i.n22() + 1) == f.n22() && i.n11() == f.n11();
  }

  inline bool E13_12(const u3::gelfand& i, const u3::gelfand& f)
  {
    return (i.n12() + 1) == f.n12() && i.n22() == f.n22() &&
           (i.n11() + 1) == f.n11();
  }

  inline bool E13_22(const u3::gelfand& i, const u3::gelfand& f)
  {
    return i.n12() == f.n12() && (i.n22() + 1) == f.n22() &&
           (i.n11() + 1) == f.n11();
  }

  // find all Gelfand states associated with an irrep
  inline std::vector<u3::gelfand> gelfandBasis(const u3::irrep& x)
  {
    std::vector<u3::gelfand> gArr;

    for (int n12 = x.n13(); n12 >= x.n23(); n12--)
      for (int n22 = x.n23(); n22 >= x.n33(); n22--)
        for (int n11 = n12; n11 >= n22; n11--)
          gArr.emplace_back(x.n13(), x.n23(), x.n33(), n12, n22, n11);

    return gArr;
  }

  // find the location (index) of a Gelfand state within the Gelfand basis above
  inline std::size_t findGelfand(const u3::gelfand& g,
                                 const std::vector<u3::gelfand>& gArr)
  {
    std::size_t i;
    for (i = 0; i < gArr.size(); i++)
      if (u3::E11(g, gArr[i])) break;
        
    return i;
  }

  /*------------------------------------------------------
  The highest weight physical basis has the following form
  | n13,n23,0  \
  | kappa, L /
  */

  /*
  // IntP(x) function, see Eq. (35) of Pan et al., Nucl. Phys. A (2016)
  inline int IntP(double x)
  {
    if (x >= 0)
      return (int)std::floor(x);
    else
      return 0;
  }

  // Multiplicity of L
  inline int Kappa(const u3::irrep& u3, int l)
  {
    return IntP((u3.n13() - u3.n33() - l + 0.0) / 2) -
           IntP((u3.n13() - u3.n23() + 1.0 - l) / 2) -
           IntP((u3.n23() - u3.n33() + 1.0 - l) / 2) + 1;
  }
  */

  inline int Kappa(const u3::irrep& ir, int L) { return su3::kmax(ir, L); }

  // alternative called 'kmax':
  inline int kmax(const u3::irrep& ir, int L) { return su3::kmax(ir, L); }

  // alternative called 'mult':
  inline int mult(const u3::irrep& ir, int L) { return su3::kmax(ir, L); }

  class physicalvector
  {
  private:
    u3::irrep irr_;
    u3::uqn_t L_;
    u3::mult_t kappa_;

  public:
    // Default constructor
    physicalvector() : kappa_(0), L_(0) { }

    // initialization
    physicalvector(const u3::irrep& irr, qcx_utils::UInt8 L) 
      : irr_(irr), 
        kappa_(qcx_utils::UInt8(Kappa(irr, L))), 
        L_(L)
    { }

    // assignment operator
    physicalvector& operator=(const physicalvector& physvec)
    {
      irr_ = physvec.irr_;
      kappa_ = physvec.kappa_;
      L_ = physvec.L_;
      return *this;
    }

    const auto& irr() const { return irr_; }

    // get SU(3)
    int lam() const { return irr_.lam(); }
    int lambda() const { return irr_.lam(); }
    int mu() const { return irr_.mu(); }

    // get U(3)
    int n13() const { return irr_.n13(); }
    int n23() const { return irr_.n23(); }
    int n33() const { return irr_.n33(); }

    // get SO(3)
    int kappa() const { return static_cast<int>(kappa_); }
    int kappamax() const { return kappa(); }
    int L() const { return static_cast<int>(L_); }
  };

  // print out a physical basis vector
  inline std::ostream& operator<<(std::ostream& s, physicalvector& physvec)
  {
    return s << "(" << physvec.lam() << "," << physvec.mu()
             << "), kappa=" << physvec.kappa() << ", L=" << physvec.L();
  }

  // find all physical states associated with an irrep
  inline std::vector<u3::physicalvector> physicalBasis(const u3::irrep& u3)
  {
    std::vector<u3::physicalvector> basis;
    for (int l = u3.n13() - u3.n33(); l >= 0; l--)
    {
      u3::physicalvector physvec(u3, l);
      if (physvec.kappa() == 0)
        continue;
      else
        basis.push_back(physvec);
    }

    return basis;
  }

} // namespace u3

#endif
