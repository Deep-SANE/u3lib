#ifndef UTILS_READWRITE_H
#define UTILS_READWRITE_H

#include <cmath>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "sparse_utils.h" // check_symmetric

namespace utils
{
  // helper function for binary writing a single value to output stream:
  template <typename T>
  inline void writeBinary(std::ostream& s, const T& t)
  {
    const char* ptr = reinterpret_cast<const char*>(&t);
    s.write(ptr, sizeof(t));
  }

  // helper function for binary writing a vector to output stream:
  template <typename T>
  inline void writeBinary(std::ostream& s, const std::vector<T>& v)
  {
    if (v.size() > 0)
      s.write(reinterpret_cast<const char*>(v.data()), v.size() * sizeof(T));
  }

  // helper function for binary reading a single value from input stream:
  template <typename T>
  inline void readBinary(std::istream& s, T& t)
  {
    char* ptr = reinterpret_cast<char*>(&t);
    s.read(ptr, sizeof(t));
  }

  // helper function for binary reading a vector from input stream:
  template <typename T>
  inline void readBinary(std::istream& s, std::vector<T>& v)
  {
    assert(v.size() > 0);
    s.read(reinterpret_cast<char*>(v.data()), v.size() * sizeof(T));
  }

  /*
  // helper function for binary writing a single value to output stream:
  template <typename T, typename V>
  inline void write(std::ostream& s, const V& v)
  {
    T t = v; // 'T' is the written type
    const char* ptr = reinterpret_cast<const char*>(&t);
    s.write(ptr, sizeof(t));
  }

  // helper function for binary reading a single value from input stream:
  template <typename T, typename V>
  inline void read(std::istream& s, V& v)
  {
    T t; // 'T' is the read type
    char* ptr = reinterpret_cast<char*>(&t);
    s.read(ptr, sizeof(t));
    v = t;
  }
  */


  // write Eigen::SparseMatrix<double> to binary file
  inline void writeBinary(const std::string& filename,
                   const Eigen::SparseMatrix<double>& mat,
                   double threshold = 1e-12)
  {
    std::ofstream out(filename, std::ios::binary);
    if (!out)
      throw std::runtime_error("Cannot open file for writing: " + filename);
      
    std::cout << "Start writing matrix to binary file" << std::endl;

    Eigen::Index rows = mat.rows();
    Eigen::Index cols = mat.cols();

    writeBinary(out, rows);
    writeBinary(out, cols);

    if (rows == cols) // symmetric case
    {
      assert(check_symmetric(mat));

      // Count number of nonzero elements above threshold
      size_t nnz = 0;
      for (int k = 0; k < mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
          if (it.row() >= it.col() && std::fabs(it.value()) > threshold)
            ++nnz;

      writeBinary(out, nnz);

      // Write entries
      for (int k = 0; k < mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
        {
          if (std::fabs(it.value()) <= threshold) continue;

          Eigen::Index i = it.row();
          Eigen::Index j = it.col();
          if (i < j) continue;
          
          writeBinary(out, i);
          writeBinary(out, j);
          writeBinary(out, it.value());
        }
    }
    else
    {
      // Count number of nonzero elements above threshold
      size_t nnz = 0;
      for (int k = 0; k < mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
          if (std::fabs(it.value()) > threshold)
            ++nnz;

      writeBinary(out, nnz);

      // Write entries
      for (int k = 0; k < mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<double>::InnerIterator it(mat, k); it; ++it)
        {
          if (std::fabs(it.value()) <= threshold) continue;

          Eigen::Index i = it.row();
          Eigen::Index j = it.col();
          
          writeBinary(out, i);
          writeBinary(out, j);
          writeBinary(out, it.value());
        }
    }
    
    std::cout << "Finish writing matrix to binary file" << std::endl;
    return;
  }

  // read Eigen::SparseMatrix<double> from binary file
  inline Eigen::SparseMatrix<double> readBinary(const std::string& filename)
  {
    std::ifstream in(filename, std::ios::binary);
    if (!in)
      throw std::runtime_error("Cannot open file for reading: " + filename);
    
    std::cout << "Start reading matrix from binary file" << std::endl;
      
    Eigen::Index rows, cols;
    size_t nnz;
    readBinary(in, rows);
    readBinary(in, cols);
    readBinary(in, nnz);

    std::vector<Eigen::Triplet<double>> triplets;

    if (rows == cols)
    {
      triplets.reserve(2 * nnz);

      for (size_t k = 0; k < nnz; ++k)
      {
        Eigen::Index i, j;
        double val;
        readBinary(in, i);
        readBinary(in, j);
        readBinary(in, val);

        triplets.emplace_back(i, j, val);
        if (i > j) triplets.emplace_back(j, i, val);
      }
    }
    else
    {
      triplets.reserve(nnz);

      for (size_t k = 0; k < nnz; ++k)
      {
        Eigen::Index i, j;
        double val;
        readBinary(in, i);
        readBinary(in, j);
        readBinary(in, val);
        triplets.emplace_back(i, j, val);
      }
    }

    Eigen::SparseMatrix<double> mat(rows, cols);
    mat.setFromTriplets(triplets.begin(), triplets.end());
      
    std::cout << "Finish reading matrix from binary file" << std::endl;
    return mat;
  }

  // write eigenvalues and eigenvectors to binary file
  inline void writeEigenBinary(const std::string& filename,
                          const std::vector<double>& evals,
                          const Eigen::MatrixXd& evecs)
  {
    std::ofstream out(filename, std::ios::binary);
    if (!out)
      throw std::runtime_error("Cannot open file for writing: " + filename);
      
    std::cout << "Start writing eigen solutions to binary file" << std::endl;

    size_t rows = evecs.rows();
    size_t cols = evecs.cols();

    writeBinary(out, rows);
    writeBinary(out, cols);

    writeBinary(out, evals);

    for (size_t i = 0; i < rows; i++)
      for (size_t j = 0; j < cols; j++)
        writeBinary(out, evecs(i, j));

    std::cout << "Finish writing eigen solutions to binary file" << std::endl;
    return;
  }

  // write eigenvalues and eigenvectors to binary file
  inline auto readEigenBinary(const std::string& filename)
  {
    std::ifstream in(filename, std::ios::binary);
    if (!in)
      throw std::runtime_error("Cannot open file for reading: " + filename);
      
    std::cout << "Start reading eigen solutions from binary file" << std::endl;

    size_t rows, cols;
    readBinary(in, rows);
    readBinary(in, cols);

    std::vector<double> evals(cols);
    readBinary(in, evals);

    Eigen::MatrixXd evecs(rows, cols);
    for (size_t i = 0; i < rows; i++)
      for (size_t j = 0; j < cols; j++)
      {
        double x;
        readBinary(in, x);
        evecs(i, j) = x;
      }

    std::cout << "Finish reading eigen solutions from binary file" << std::endl;
    return std::make_pair(evals, evecs);
  }

} // namespace utils

#endif
