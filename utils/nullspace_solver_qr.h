#ifndef UTILS_NULLSPACE_SOLVER_QR_H
#define UTILS_NULLSPACE_SOLVER_QR_H

#define FIX_MULTIS_AMBIGUITY

#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "nullspace_error.h"
#include "sparse_utils.h"

namespace utils
{
  // Solving an upper square matrix by backward substitution
  inline Eigen::MatrixXd BackwardSubstitution(const Eigen::MatrixXd& mat,
                                              int nullity)
  {
    int row_num = mat.cols();
    double epsilon = 1e-10;

    Eigen::MatrixXd result(row_num, nullity);

    for (int i = 0; i < nullity; i++)
    {
      for (int j = row_num - i; j < row_num; j++)
        result(j, i) = 0;

      result(row_num - i - 1, i) = 1.0;

      for (int j = row_num - i - 2; j >= 0; j--)
      {
        double sum = 0.0;

        for (int k = j + 1; k < row_num; k++)
          sum += mat(j, k) * result(k, i);

        if (std::fabs(mat(j, j)) < epsilon)
          result(j, i) = 0.0;

        else
          result(j, i) = -sum / mat(j, j);
      }
    }

    return result;
  }

  // A function to perform the modified Gram-Schmidt orthonormalization -- a
  // more stable version of the traditional algorithm
  inline void GramSchmidt(Eigen::MatrixXd& mat, double tol=1e-12)
  {
    int col_num = mat.cols();

    for (int i = 0; i < col_num; i++)
    {
      // Orthogonalize against previous columns
      for (int j = 0; j < i; j++)
      {
        double proj = mat.col(j).dot(mat.col(i));
        mat.col(i) -= proj * mat.col(j);
      }

      // Normalize
      double n = mat.col(i).norm();

      if (n > tol) mat.col(i) /= n;
      else mat.col(i).setZero();
    }

    // traditional version:
    /*
    matMGS = mat;

    for (int j = 0; j < matMGS.cols(); j++)
    {
      // Orthogonalization
      for (int k = 0; k < j; k++)
      {
        double projection = 0.0;
        for (int i = 0; i < matMGS.rows(); i++)
          projection += matMGS(i, j) * matMGS(i, k);
        for (int i = 0; i < matMGS.rows(); i++)
          matMGS(i, j) -= projection * matMGS(i, k);
      }
      matMGS.col(j) /= matMGS.col(j).norm();
    }

    return matMGS;
    */
  }

  // Reduced row echelon form to bring row-equivalent matricies to the canonical form
  inline void ReducedRowEchelon(Eigen::MatrixXd& mat, double tol=1e-12)
  {
    int rows = mat.rows();
    int cols = mat.cols();

    int r = 0; // pivot row

    for (int c = 0; c < cols && r < rows; c++)
    {
      // Find pivot
      int pivot = r;
      double maxval = std::abs(mat(r,c));

      for (int i = r + 1; i < rows; i++)
      {
        double val = std::abs(mat(i, c));
        if (val > maxval)
        {
          maxval = val;
          pivot = i;
        }
      }

      // Skip if pivot ~ zero
      if (maxval < tol) continue;

      // Swap rows
      if (pivot != r) mat.row(r).swap(mat.row(pivot));

      // Normalize pivot row
      double piv = mat(r, c);
      mat.row(r) /= piv;

      // Eliminate all other rows
      for (int i = 0; i < rows; i++)
      {
        if (i == r) continue;

        double factor = mat(i, c);

        if (std::fabs(factor) > tol)
        {
          mat.row(i) -= factor * mat.row(r);
          mat(i,c) = 0.0;
        }
      }

      r++;
    }

    // Clean small numbers
    mat = (mat.array().abs() < tol).select(0.0, mat);
  }

  inline void SignConvention(Eigen::MatrixXd& mat, double tol=1e-15)
  {
    for (int i = 0; i < mat.cols(); i++)
    {
      int j;
      for (j = 0; j < mat.rows(); j++)
        if (std::fabs(mat(j, i)) > tol) break;

      if (mat(j, i) < 0) mat.col(i) *= (-1.0);

      for (j = 0; j < mat.rows(); j++)
        if (std::fabs(mat(j, i)) <= tol) mat(j, i) = 0.0;
    }
  }

  // Solving the null space of a rectangular matrix
  /*
  The main strategy to solve the null space of a matrix A of size m x n is as
  follows i) Apply QR decomposition on the matrix: AP = QR, ii) Matrix R has
  size m x n, thus it will be resized to n x n, iii) Solve matrix R by backward
  substitution to obtain the null space of A, iv) However, columns of A has been
  permuted, the last step is to permute the components of the null space vectors
  back by multiplying the permutation matrix P with the null space matrix.
  */
  inline Eigen::MatrixXd null_space(const Eigen::MatrixXd& mat)
  {
    if (mat.size() <= 1)
    {
      Eigen::Matrix<double, 1, 1> one;
      one(0, 0) = 1;
      return one;
    }
    else
    {
      // Create a matrix to store the null space
      Eigen::MatrixXd nullSpace;

      // Sparse matrix solver
      if (sparsity(mat) >= 0.7)
      {
        // Convert Matrix to SparseMatrix
        Eigen::SparseMatrix<double> spar_mat = mat.sparseView();
        spar_mat.makeCompressed();

        // Create a SparseQR decomposition of the projection matrix
        Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
            qr;

        // Compute the qR decomposition of the projection matrix
        qr.setPivotThreshold(1e-6);
        qr.compute(spar_mat);

        // Get the upper triangular matrix and resize it to square
        Eigen::MatrixXd R(qr.matrixR());

        assert((R.rows() == mat.rows()) && (R.cols() == mat.cols()));

        if (mat.cols() > mat.rows())
        {
          R.conservativeResize(mat.cols(), mat.cols());
          R.block(mat.rows(), 0, mat.cols() - mat.rows(), mat.cols()).setZero();
        }

        // Get the nullity -- dimension of the nullspace
        int nullity = mat.cols() - qr.rank();

        nullSpace = BackwardSubstitution(R, nullity);
        nullSpace = qr.colsPermutation() * nullSpace;
      }
      else
      {
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;

        // Compute the qR decomposition of the projection matrix
        // qr.setThreshold(1e-16);
        qr.compute(mat);

        // Get the upper triangular matrix and resize it to square
        Eigen::MatrixXd R(qr.matrixR().template triangularView<Eigen::Upper>());

        assert((R.rows() == mat.rows()) && (R.cols() == mat.cols()));

        if (mat.cols() > mat.rows())
        {
          R.conservativeResize(mat.cols(), mat.cols());
          R.block(mat.rows(), 0, mat.cols() - mat.rows(), mat.cols()).setZero();
        }

        // Get the nullity -- dimension of the nullspace
        int nullity = qr.dimensionOfKernel();

        nullSpace = BackwardSubstitution(R, nullity);
        nullSpace = qr.colsPermutation() * nullSpace;
      }

#ifdef FIX_MULTIS_AMBIGUITY
      Eigen::MatrixXd nullSpace_t = nullSpace.transpose();
      ReducedRowEchelon(nullSpace_t);
      nullSpace = nullSpace_t.transpose();
#endif

      GramSchmidt(nullSpace);
      SignConvention(nullSpace);

      assert(nullspace_error(mat, nullSpace) < 1e-6);

      return nullSpace;
    }
  }

  // "direct" version for sparse input matrices
  inline Eigen::MatrixXd null_space(const Eigen::SparseMatrix<double>& mat)
  {
    if (mat.size() <= 1)
    {
      Eigen::Matrix<double, 1, 1> one;
      one(0, 0) = 1;
      return one;
    }
    else
    {
      // Create a SparseQR decomposition of the projection matrix
      Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
          qr;

      // Compute the qR decomposition of the projection matrix
      qr.setPivotThreshold(1e-6);
      qr.compute(mat);

      // Get the upper triangular matrix and resize it to square
      Eigen::MatrixXd R(qr.matrixR());

      assert((R.rows() == mat.rows()) && (R.cols() == mat.cols()));

      if (mat.cols() > mat.rows())
      {
        R.conservativeResize(mat.cols(), mat.cols());
        R.block(mat.rows(), 0, mat.cols() - mat.rows(), mat.cols()).setZero();
      }

      // Get the nullity -- dimension of the nullspace
      int nullity = mat.cols() - qr.rank();

      Eigen::MatrixXd res = BackwardSubstitution(R, nullity);
      res = qr.colsPermutation() * res;

#ifdef FIX_MULTIS_AMBIGUITY
      Eigen::MatrixXd res_t = res.transpose();
      ReducedRowEchelon(res_t);
      res = res_t.transpose();
#endif

      GramSchmidt(res);
      SignConvention(res);
      
      assert(nullspace_error(mat, res) < 1e-8);

      return res;
    }
  }
} // namespace utils

#endif
