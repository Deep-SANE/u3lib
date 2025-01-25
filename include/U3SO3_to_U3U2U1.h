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

  /*------------------------------------------------------
  The gelfand (canonical) basis states for each M >= 0 can be expressed by two
  parameters q,t
  |  n13,n23,0 \
  |  M+2*q-t,t  >
  |     M+q    /
  where
  0 <= M <= n13
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

  std::vector<int> expand_to_gelfand(int M, std::vector<u3::gelfand>& gArr)
  {
    std::vector<int> indices;

    for (int i = 0; i < gArr.size(); i++)
    {
      if (M == gArr[i].e11() - gArr[i].e22())
        indices.push_back(i);
      /*
      int q = gArr[i].n11 - M;
      int qmin = std::max(
          gArr[i].n22,
          IntM((gArr[i].n22 - M + gArr[i].n23) / 2));
      int qmax = (int)std::floor(
          (gArr[i].n13 - M + gArr[i].n22) / 2);

      if ((M + 2 * q) == (gArr[i].n12 + gArr[i].n22) && q <= qmax &&
          q >= qmin)
        indices.push_back(i);
      */
    }
    return indices;
  }

  /*------------------------------------------------------
  Finding the transformation by solving a system of linear equations resulted
  from acting L+ = sqrt(2)*(E13+E32)
  */
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  Lprojection(u3::physicalvector& physvec, std::vector<u3::gelfand>& gArr,
              std::vector<int>& indices)
  {

    std::vector<int> new_indices;

    for (int i = 0; i < gArr.size(); i++)
    {

      for (int j = 0; j < indices.size(); j++)
      {
        // action of E13
        if (u3::E13_12(gArr[indices[j]], gArr[i]))
        {
          new_indices.push_back(i);
          break;
        }
        else if (u3::E13_22(gArr[indices[j]], gArr[i]))
        {
          new_indices.push_back(i);
          break;
        }

        // action of E32
        else if (u3::E23_12(gArr[i], gArr[indices[j]]))
        {
          new_indices.push_back(i);
          break;
        }
        else if (u3::E23_22(gArr[i], gArr[indices[j]]))
        {
          new_indices.push_back(i);
          break;
        }
      }
    }

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> projector(
        new_indices.size(), indices.size());

    for (int i = 0; i < projector.rows(); i++)
    {
      for (int j = 0; j < projector.cols(); j++)
      {

        // action of E13
        if (u3::E13_12(gArr[indices[j]], gArr[new_indices[i]]))
          projector(i, j) = gArr[indices[j]].e13_12();
        else if (u3::E13_22(gArr[indices[j]], gArr[new_indices[i]]))
          projector(i, j) = gArr[indices[j]].e13_22();

        // action of E32
        else if (u3::E23_12(gArr[new_indices[i]], gArr[indices[j]]))
          projector(i, j) = gArr[indices[j]].e32_12();
        else if (u3::E23_22(gArr[new_indices[i]], gArr[indices[j]]))
          projector(i, j) = gArr[indices[j]].e32_22();

        else
          projector(i, j) = 0.0;
      }
    }
    return projector;
  }

  // get the transformation
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  phys_to_cano(u3::physicalvector& physvec, std::vector<int>& indices,
               std::vector<u3::gelfand>& gArr)
  {
    indices = u3::expand_to_gelfand(physvec.L, gArr);

    // Building the projection matrix
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Lproj =
        u3::Lprojection(physvec, gArr, indices);
    // and solving its null space
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nullSpace =
        null_space(Lproj);
    // after solving the nullspace, the dimension of the nullspace is the
    // multiplicity of L, and the multiplicity within the physical basis state
    // labeling is modified

    // return the multipplicity
    physvec.kappa = nullSpace.cols();

    return nullSpace;
  }

} // namespace u3

#endif

/*
// Building the projection matrix from the physical basis
// |[n13,n23,n33],kappa,L> onto gelfand basis

// The original version
std::vector<u3::gelfand> expand_to_gelfand(physicalvector& physvec)
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
      if (u3::E11(g, newgArr[j]))
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
      if (u3::E13_12(gArr[j],newgArr[i]))
      {
        projector(i, j) = gArr[j].e13_12();
      }
      else if (u3::E13_22(gArr[j],newgArr[i]))
      {
        projector(i, j) = gArr[j].e13_22();
      }

      // action of E32
      else if (u3::E23_12(newgArr[i],gArr[j]))
      {
        projector(i, j) = gArr[j].e32_12();
      }
      else if (u3::E23_22(newgArr[i],gArr[j]))
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



// The version that builds the angular momentum projection matrix directly
Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
Lprojection2(physicalvector& physvec, std::vector<u3::gelfand>& gArr,
std::vector<int>& indices)
{

  int row_num=1;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
projector(indices.size(), indices.size());

  for (int i = 0; i < gArr.size(); i++)
  {
    for (int j=0; j<indices.size(); j++)
    {
      // action of E13
      if (u3::E13_12(gArr[indices[j]],gArr[i]))
      {
        projector(row_num-1, j) = gArr[indices[j]].e13_12();
        for (int k=j+1; k<indices.size(); k++)
        {
          if (u3::E13_22(gArr[indices[k]],gArr[i]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e13_22();
          }
          else if (u3::E23_12(gArr[i],gArr[indices[k]]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e32_12();
          }
          else if (u3::E23_22(gArr[i],gArr[indices[k]]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e32_22();
          }
        }
        row_num++;
        break;
      }
      else if (u3::E13_22(gArr[indices[j]],gArr[i]))
      {
        projector(row_num-1, j) = gArr[indices[j]].e13_22();
        for (int k=j+1; k<indices.size(); k++)
        {
          if (u3::E13_12(gArr[indices[k]],gArr[i]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e13_12();
          }
          else if (u3::E23_12(gArr[i],gArr[indices[k]]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e32_12();
          }
          else if (u3::E23_22(gArr[i],gArr[indices[k]]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e32_22();
          }
        }
        row_num++;
        break;
      }

      // action of E32
      else if (u3::E23_12(gArr[i],gArr[indices[j]]))
      {
        projector(row_num-1, j) = gArr[indices[j]].e32_12();
        for (int k=j+1; k<indices.size(); k++)
        {
          if (u3::E13_12(gArr[indices[k]],gArr[i]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e13_12();
          }
          else if (u3::E13_22(gArr[indices[k]],gArr[i]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e13_22();
          }
          else if (u3::E23_22(gArr[i],gArr[indices[k]]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e32_22();
          }
        }
        row_num++;
        break;
      }
      else if (u3::E23_22(gArr[i],gArr[indices[j]]))
      {
        projector(row_num-1, j) = gArr[indices[j]].e32_22();
        for (int k=j+1; k<indices.size(); k++)
        {
          if (u3::E13_12(gArr[indices[k]],gArr[i]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e13_12();
          }
          else if (u3::E13_22(gArr[indices[k]],gArr[i]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e13_22();
          }
          else if (u3::E23_12(gArr[i],gArr[indices[k]]))
          {
            projector(row_num-1, k) = gArr[indices[k]].e32_12();
          }
        }
        row_num++;
        break;
      }
      else
      {
        projector(row_num-1, j) = 0.0;
      }
    }
  }

  projector.conservativeResize(row_num-1, indices.size());
  return projector;
}
*/
