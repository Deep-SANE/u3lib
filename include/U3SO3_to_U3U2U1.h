// protection against multiple header file inclusion
#ifndef U3SO3_to_U3U2U1_H
#define U3SO3_to_U3U2U1_H

#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>

#include "U3basis.h"
#include "utils/nullspaceSolver.h"

namespace u3
{

  // construction of the physical basis: the physical basis vectors will be
  // stored in an array
  std::vector<physicalvector> physicalbasis(u3::irrep& u3)
  {
    std::vector<physicalvector> basis;
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

  /*------------------------------------------------------
  The gelfand (canonical) basis state for each L can be expressed by two
  parameters q,t
  |  n13,n23,0 \
  |  L+2*q-t,t  >
  |     L+q    /
  where
  0 <= L <= n13
  0 <= t <= n23
  max(t,round((t-L+n23)/2)) <= q <= Int((n13-L+t)/2)
  */

  // Definition of IntM(x) funtion
  int IntM(double x)
  {
    if ((x - std::floor(x)) == 0.0)
      return (int)std::floor(x);
    else
      return (int)std::floor(x) + 1;
  }

  // construction of the canonical basis: the canonical basis vectors will be
  // stored in an array
  std::vector<u3::gelfand> canonicalbasis(physicalvector& physvec)
  {
    std::vector<u3::gelfand> gArr;
    for (int t = (physvec.toprow.n23 - physvec.toprow.n33); t >= 0; t--)
    {
      for (int q = (int)std::floor(
               (physvec.toprow.n13 - physvec.toprow.n33 - physvec.L + t) / 2);
           q >= std::max(t, IntM((t - physvec.L + physvec.toprow.n23 -
                                  physvec.toprow.n33) /
                                 2));
           q--)
      {
        u3::gelfand g(physvec.toprow.n13, physvec.toprow.n23,
                      physvec.toprow.n33, physvec.L + 2 * q - t, t,
                      physvec.L + q);
        if (g.check())
          gArr.push_back(g);
      }
    }
    return gArr;
  }

  /*------------------------------------------------------
  Finding the transformation by solving a system of linear equations resulted
  from acting L+ = sqrt(2)*(E13+E32)
  */
  void newstate(u3::gelfand& g, std::vector<u3::gelfand>& newgArr)
  {
    if (newgArr.size() == 0)
    {
      newgArr.push_back(g);
    }
    else
    {
      int flag = 0;
      for (int j = 0; j < newgArr.size(); j++)
      {
        if (u3::compare(g, newgArr[j]))
        {
          flag = 1;
          break;
        }
      }
      if (flag == 0)
        newgArr.push_back(g);
    }
    return;
  }

  // Building the projection matrix from the physical basis
  // |[n13,n23,n33],kappa,L> onto gelfand basis
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  Lprojection(physicalvector& physvec, std::vector<u3::gelfand>& gArr)
  {

    std::vector<u3::gelfand> newgArr;

    for (int i = 0; i < gArr.size(); i++)
    {

      // determine the action of E13
      u3::gelfand g1 = gArr[i];
      g1.n12++;
      g1.n11++;
      if (g1.check())
        newstate(g1, newgArr);

      u3::gelfand g2 = gArr[i];
      g2.n22++;
      g2.n11++;
      if (g2.check())
        newstate(g2, newgArr);

      // determine the action of E32
      u3::gelfand g3 = gArr[i];
      g3.n12--;
      if (g3.check())
        newstate(g3, newgArr);
      u3::gelfand g4 = gArr[i];
      g4.n22--;
      if (g4.check())
        newstate(g4, newgArr);
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> projector(
        newgArr.size(), gArr.size());

    for (int i = 0; i < projector.rows(); i++)
    {
      for (int j = 0; j < projector.cols(); j++)
      {

        // action of E13
        if (newgArr[i].n12 == (gArr[j].n12 + 1) &&
            newgArr[i].n22 == gArr[j].n22 &&
            newgArr[i].n11 == (gArr[j].n11 + 1))
        {
          projector(i, j) = gArr[j].e13_12();
        }
        else if (newgArr[i].n12 == gArr[j].n12 &&
                 newgArr[i].n22 == (gArr[j].n22 + 1) &&
                 newgArr[i].n11 == (gArr[j].n11 + 1))
        {
          projector(i, j) = gArr[j].e13_22();
        }

        // action of E32
        else if (newgArr[i].n12 == (gArr[j].n12 - 1) &&
                 newgArr[i].n22 == gArr[j].n22 && newgArr[i].n11 == gArr[j].n11)
        {
          projector(i, j) = gArr[j].e32_12();
        }
        else if (newgArr[i].n12 == gArr[j].n12 &&
                 newgArr[i].n22 == (gArr[j].n22 - 1) &&
                 newgArr[i].n11 == gArr[j].n11)
        {
          projector(i, j) = gArr[j].e32_22();
        }

        else
        {
          projector(i, j) = 0.0;
        }
      }
    }
    return projector;
  }

  // Solving the null space of the projection matrix
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  phys_to_cano(physicalvector& physvec,
               Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& Lproj)
  {

    // Create a matrix to store the null space
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nullSpace =
        null_space(Lproj);

    // return the multipplicity
    physvec.kappa = nullSpace.cols();

    return nullSpace;
  }

} // namespace u3

#endif
