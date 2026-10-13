#ifndef UTILS_SPARSE_UTILS_H
#define UTILS_SPARSE_UTILS_H

#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace utils
{
  // Determine sparsity of a matrix
  inline double sparsity(const Eigen::MatrixXd& mat)
  {
    double s = 0;
    for (int i = 0; i < mat.rows(); i++)
    {
      for (int j = 0; j < mat.cols(); j++)
        if (mat(i, j) == 0.0)
          s++;
    }
    return s / (mat.rows() * mat.cols());
  }

  // Check whether a sparse matrix is (numerically) symmetric.
  inline bool check_symmetric(const Eigen::SparseMatrix<double>& SA)
  {
    Eigen::MatrixXd A = SA;

    const auto m = A.rows();
    const auto n = A.cols();

    for (size_t i = 0; i < n; i++)
      for (size_t j = 0; j < i; j++)
        if (std::fabs(A(i, j) - A(j, i)) > 1e-4)
          return false;

    return true;
  }
} // namespace utils

#endif
