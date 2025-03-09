#ifndef U3LIB_U3SO3_to_U3U2U1_H
#define U3LIB_U3SO3_to_U3U2U1_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "U3basis.h"
#include "utils/nullspaceSolver.h"

namespace u3
{

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

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  lowerM(int kappa_max, int L, int M, std::vector<u3::gelfand>& gArr,
         Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& transformation,
         std::vector<int>& canoBasis, std::vector<int>& canoBasis_new)
  {
    // normalize coefficient
    double n = std::sqrt(2.0 / ((L - M) * (L + M + 1.0)));
    // input M is for the new state

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> transformation_new(
        canoBasis_new.size(), kappa_max);
    transformation_new.setZero();

    // basis states are obtained by acting with two generators E31 and E23
    for (int i = 0; i < canoBasis_new.size(); i++)
    {
      for (int j = 0; j < canoBasis.size(); j++)
      {
        if (u3::E13_12(gArr[canoBasis_new[i]], gArr[canoBasis[j]]))
          for (int k = 0; k < kappa_max; k++)
            transformation_new(i, k) +=
                n * transformation(j, k) * gArr[canoBasis[j]].e31_12();
        else if (u3::E13_22(gArr[canoBasis_new[i]], gArr[canoBasis[j]]))
          for (int k = 0; k < kappa_max; k++)
            transformation_new(i, k) +=
                n * transformation(j, k) * gArr[canoBasis[j]].e31_22();
        else if (u3::E23_12(gArr[canoBasis[j]], gArr[canoBasis_new[i]]))
          for (int k = 0; k < kappa_max; k++)
            transformation_new(i, k) +=
                n * transformation(j, k) * gArr[canoBasis[j]].e23_12();
        else if (u3::E23_22(gArr[canoBasis[j]], gArr[canoBasis_new[i]]))
          for (int k = 0; k < kappa_max; k++)
            transformation_new(i, k) +=
                n * transformation(j, k) * gArr[canoBasis[j]].e23_22();
      }

      for (int k = 0; k < kappa_max; k++)
        if (std::fabs(transformation_new(i, k)) < 1e-12)
          transformation_new(i, k) = 0.0;
    }
    return transformation_new;
  }

  // get the transformation for M <= L
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  phys_to_cano(u3::physicalvector& physvec, int M, std::vector<int>& canoBasis,
               std::vector<u3::gelfand>& gArr)
  {
    canoBasis = u3::expand_to_gelfand(physvec.L, gArr);

    // Building the projection matrix
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Lproj =
        u3::Lprojection(physvec, gArr, canoBasis);
    // and solving its null space
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> transformation =
        null_space(Lproj);

    // ensure the nullspace is solved correctly
    assert(physvec.kappa == transformation.cols());

    for (int m = physvec.L - 1; m >= std::fabs(M); m--)
    {
      // store new expansion
      std::vector<int> canoBasis_new = u3::expand_to_gelfand(m, gArr);

      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> transformation_new =
          lowerM(physvec.kappa, physvec.L, m, gArr, transformation, canoBasis,
                 canoBasis_new);

      transformation = transformation_new;
      canoBasis = canoBasis_new;
    }

    if (M < 0)
    {
      canoBasis = u3::expand_to_gelfand(M, gArr);

      for (int j = 0; j < canoBasis.size(); j++)
      {
        int phase =
            ((physvec.toprow.lam() + physvec.L + gArr[canoBasis[j]].n22()) %
                 2 ==
             0)
                ? 1
                : -1;
        transformation.row(j) *= phase;
      }
    }

    return transformation;
  }

} // namespace u3

#endif
