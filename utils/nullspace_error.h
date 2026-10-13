#ifndef UTILS_NULLSPACE_ERROR_H
#define UTILS_NULLSPACE_ERROR_H

#include <cassert>
#include <cmath>

#include <Eigen/Dense>

namespace utils
{
  // A function to return the standard deviation (error) of solving the
  // nullspace, it is calculated by going through all elements of
  // (projection_matrix * nullspace) and determine the standard deviation
  //
  // T can be either a dense or a sparse Eigen matrix type
  template <typename T>
  double nullspace_error(const T& mat, const Eigen::MatrixXd& nullspace)
  {
    assert(mat.cols() == nullspace.rows());

    // nothing to check if the null space is empty
    if (nullspace.cols() == 0)
      return 0.0;

    /*
    if (mat.size() <= 1)
      return 0.0;
    */

    assert(mat.cols() == nullspace.rows());

    Eigen::MatrixXd zeros = mat * nullspace;

    double error = 0.0;

    /*
    for (int i = 0; i < zeros.rows(); i++)
      for (int j = 0; j < zeros.cols(); j++)
        error += std::fabs(zeros(i, j));
    error /= zeros.size();
    */

    // error = maximum of norms of residual vectors:
    for (int i = 0; i < zeros.cols(); i++)
    {
      auto norm = zeros.col(i).norm();
      if (norm > error)
        error = norm;
    }
    // make error relative their lentgh:
    error /= std::sqrt(zeros.cols());

    return error;
  }
} // namespace utils

#endif
