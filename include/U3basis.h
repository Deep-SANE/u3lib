#ifndef U3LIB_U3BASIS_H
#define U3LIB_U3BASIS_H

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include "config.h"

namespace u3
{

  /*------------------------------------------------------
  A structure for U(3) irrep [n13,n23,n33]
  */
  struct irrep
  {
    int n13, n23, n33;

    // Default constructor
    irrep() : n13(0), n23(0), n33(0) { }

    // initialization
    irrep(int n1, int n2, int n3)
    {
      n13 = n1;
      n23 = n2;
      n33 = n3;
    }

    // construction from Elliott's SU(3) (lambda,mu)
    irrep(int lambda, int mu)
    {
      n13 = lambda + mu;
      n23 = mu;
      n33 = 0;
    }

    // get lambda and mu
    int lam() { return n13 - n23; }
    int mu() { return n23 - n33; }

    // assignment operator
    irrep& operator=(const irrep& u3)
    {
      n13 = u3.n13;
      n23 = u3.n23;
      n33 = u3.n33;
      return *this;
    }

    // dimension
    int dim()
    {
      return (1 + n13 - n23) * (1 + n23 - n33) * (2 + n13 - n33) / 2;
    }

    // check order of quantum numbers
    bool check()
    {
      if (n13 >= n23 && n23 >= n33)
        return true;
      else
        return false;
    }

    // normalize the irrep
    void normalize()
    {
      n13 -= n33;
      n23 -= n33;
      n33 = 0;
    }
  };

  // print u3irrep
  inline std::ostream& operator<<(std::ostream& s, irrep& u3)
  {
    return s << "[" << u3.n13 << "," << u3.n23 << "," << u3.n33 << "]";
  }

  // print gelfand state
  inline std::ostream& operator<<(std::ostream& s, gelfand& g)
  {
    return s << g.n13() << "," << g.n23() << "," << g.n33() << ";" << g.n12()
             << "," << g.n22() << ";" << g.n11();
  }

  // compare two gelfand patterns of the same irrep
  bool E11(const u3::gelfand& i, const u3::gelfand& f)
  {
    return i.n12() == f.n12() && i.n22() == f.n22() && i.n11() == f.n11();
  }

  // define action of the raising generators, that of the lowering generators
  // is simply swapping the initial and final states
  bool E12(const u3::gelfand& i, const u3::gelfand& f)
  {
    return i.n12() == f.n12() && i.n22() == f.n22() && (i.n11() + 1) == f.n11();
  }

  bool E23_12(const u3::gelfand& i, const u3::gelfand& f)
  {
    return (i.n12() + 1) == f.n12() && i.n22() == f.n22() && i.n11() == f.n11();
  }

  bool E23_22(const u3::gelfand& i, const u3::gelfand& f)
  {
    return i.n12() == f.n12() && (i.n22() + 1) == f.n22() && i.n11() == f.n11();
  }

  bool E13_12(const u3::gelfand& i, const u3::gelfand& f)
  {
    return (i.n12() + 1) == f.n12() && i.n22() == f.n22() && (i.n11() + 1) == f.n11();
  }

  bool E13_22(const u3::gelfand& i, const u3::gelfand& f)
  {
    return i.n12() == f.n12() && (i.n22() + 1) == f.n22() && (i.n11() + 1) == f.n11();
  }

  // find all Gelfand states associated with an irrep
  std::vector<u3::gelfand> gelfandBasis(u3::irrep& x)
  {
    std::vector<u3::gelfand> gArr;

    for (int n12 = x.n13; n12 >= x.n23; n12--)
      for (int n22 = x.n23; n22 >= x.n33; n22--)
        for (int n11 = n12; n11 >= n22; n11--)
          gArr.emplace_back(x.n13, x.n23, x.n33, n12, n22, n11);

    return gArr;
  }

  // find the location (index) of a Gelfand state within the Gelfand basis above
  std::size_t findGelfand(const u3::gelfand& g, const std::vector<u3::gelfand>& gArr)
  {
    std::size_t i;

    for (i = 0; i < gArr.size(); i++)
      if (u3::E11(g, gArr[i]))
        break;

    return i;
  }

  /*------------------------------------------------------
  The highest weight physical basis has the following form
  | n13,n23,0  \
  | kappa, L /
  */

  // IntP(x) function, see Eq. (35) of Pan et al., Nucl. Phys. A (2016)
  int IntP(double x)
  {
    if (x >= 0)
      return (int)std::floor(x);
    else
      return 0;
  }

  // Multiplicity of L
  int Kappa(u3::irrep& u3, int l)
  {
    return IntP((u3.n13 - u3.n33 - l + 0.0) / 2) -
           IntP((u3.n13 - u3.n23 + 1.0 - l) / 2) -
           IntP((u3.n23 - u3.n33 + 1.0 - l) / 2) + 1;
  }

  struct physicalvector
  {
    u3::irrep toprow;
    int L;
    int kappa;

    // Default constructor
    physicalvector() : kappa(0), L(0) { }

    // initialization
    physicalvector(u3::irrep u3, int l)
    {
      toprow = u3;
      kappa = Kappa(u3, l);
      L = l;
    }

    // assignment operator
    physicalvector& operator=(const physicalvector& physvec)
    {
      toprow = physvec.toprow;
      kappa = physvec.kappa;
      L = physvec.L;
      return *this;
    }
  };

  // print out a physical basis vector
  inline std::ostream& operator<<(std::ostream& s, physicalvector& physvec)
  {
    return s << "(" << physvec.toprow.lam() << "," << physvec.toprow.mu()
             << "), kappa=" << physvec.kappa << ", L=" << physvec.L;
  }

  // find all physical states associated with an irrep
  std::vector<u3::physicalvector> physicalBasis(u3::irrep& u3)
  {
    std::vector<u3::physicalvector> basis;
    int l;

    for (l = u3.n13 - u3.n33; l >= 0; l--)
    {
      u3::physicalvector physvec(u3, l);
      if (physvec.kappa == 0)
        continue;
      else
        basis.push_back(physvec);
    }

    return basis;
  }

} // namespace u3

#endif
