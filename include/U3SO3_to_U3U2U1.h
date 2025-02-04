#ifndef U3LIB_U3SO3_to_U3U2U1_H
#define U3LIB_U3SO3_to_U3U2U1_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

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
      if (M == gArr[i].e11() - gArr[i].e22())
        indices.push_back(i);

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
