#ifndef U3LIB_U3U2XU3U2_H
#define U3LIB_U3U2XU3U2_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#ifdef U3LIB_TIME_LOG
#include <chrono>
#include <iostream>
#endif

#include <Eigen/Dense>

#ifdef U3LIB_SPARSE_MAT
#include <Eigen/Sparse>
#endif

#include "cgc_so3/CGC_SO3.h"
#include "decomposition.h"
#include "U3basis.h"
#include "utils/inverseSolver.h"
#include "utils/nullspaceSolver.h"

extern "C"
{
#include <wigxjpf.h>
}

// formula for mapping states to coupling coefficients
// there are three irreps u3a, u3b, u3c
// the coupling coefficients are stored in a matrix whose number of columns is
// the outer multiplicity rho and number of rows is dim(u3a)*dim(u3b)*dim(u3c)
// for matching the Gelfand states associated with each irrep to the row in the
// coefficient matrix, use the following formula: row = (ia*dim(u3b) +
// ib)*dim(u3c) + ic provided that the loops are nested in this order: ia -> ib
// -> ic

namespace u3
{
  class WigU3U2_coeffs
  {
    u3::irrep u3a_, u3b_, u3c_;

    std::size_t rhomax_;

    std::size_t D1_, D2_, D3_, D4_;

    struct Record
    {
      int p1, q1;
      int kmin, kmax;
      std::size_t first_row;
    };

    std::vector<Record> records_;

    std::vector<int> computed_;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        coeffs_;

    std::size_t record_index(int pa, int qa, int pb, int qb) const
    {
      assert(pa >= 0);
      assert(pa <= u3a_.lam());

      assert(qa >= 0);
      assert(qa <= u3a_.mu());

      assert(pb >= 0);
      assert(pb <= u3b_.lam());

      assert(qb >= 0);
      assert(qb <= u3b_.mu());

      std::size_t index = ((pa * D2_ + qa) * D3_ + pb) * D4_ + qb;

      assert(index < records_.size());

      return index;
    }

    const auto& record(int pa, int qa, int pb, int qb) const
    {
      return records_[record_index(pa, qa, pb, qb)];
    }

    auto& record(int pa, int qa, int pb, int qb)
    {
      return records_[record_index(pa, qa, pb, qb)];
    }

  public:
    struct Accessor
    {
      bool allowed;
      bool computed;
      std::size_t row;
    };

    WigU3U2_coeffs(const u3::irrep& u3a, const u3::irrep& u3b,
                   const u3::irrep& u3c, int rhomax)
        : u3a_(u3a), u3b_(u3b), u3c_(u3c), rhomax_(rhomax)
    {
      assert(rhomax > 0);

      D1_ = u3a.lam() + 1;
      D2_ = u3a.mu() + 1;
      D3_ = u3b.lam() + 1;
      D4_ = u3b.mu() + 1;

      records_.resize(D1_ * D2_ * D3_ * D4_);

      int temp = 2 * (u3c.n13() - 2 * u3c.n23() + u3c.n33());

      std::size_t row = 0;

      for (int pa = 0; pa <= u3a.lam(); pa++)
        for (int qa = 0; qa <= u3a.mu(); qa++)
        {
          u3::gelfand gA(u3a.n13(), u3a.n23(), u3a.n33(), u3a.n23() + pa,
                         u3a.n23() - qa, u3a.n23() - qa);

          for (int pb = 0; pb <= u3b.lam(); pb++)
            for (int qb = 0; qb <= u3b.mu(); qb++)
            {
              u3::gelfand gB(u3b.n13(), u3b.n23(), u3b.n33(), u3b.n23() + pb,
                             u3b.n23() - qb, u3b.n23() - qb);

              int Q = temp - gA.epsilon() - gB.epsilon();
              assert(Q % 3 == 0);
              int R = Q / 3;

              int p1 = (R >= 0) ? R : 0;
              int q1 = (R >= 0) ? 0 : -R;

              // epsilon-selection rule:
              int ub1 = std::min(u3c.lam() - p1, u3c.mu() - q1);

              // Lambda-selection rule:
              int absR = std::abs(R);
              int lb2 = std::abs(gA.Lambdax2() - gB.Lambdax2()) - absR;
              int ub2 = gA.Lambdax2() + gB.Lambdax2() - absR;

              // combine selection rules:
              int kmin = (std::max(0, lb2) + 1) / 2;
              int kmax = std::min(2 * ub1, ub2) / 2;

              record(pa, qa, pb, qb) = {p1, q1, kmin, kmax, row};

              if (kmax >= kmin)
                row += kmax - kmin + 1;
            }
        }

      coeffs_.resize(row, rhomax);
      coeffs_.setZero();

      computed_.resize(row, 0);
    }

    void normalize()
    {
      u3a_.normalize();
      u3b_.normalize();
      u3c_.normalize();
    }

    Accessor access(int n12a, int n22a, int n12b, int n22b, int n12c,
                    int n22c) const
    {
      int pa = n12a - u3a_.n23();
      int qa = u3a_.n23() - n22a;
      int pb = n12b - u3b_.n23();
      int qb = u3b_.n23() - n22b;

      const auto& rec = record(pa, qa, pb, qb);

      if (rec.kmax < rec.kmin)
        return {false, false, 0};

      int k1 = n12c - u3c_.n23() - rec.p1;
      int k2 = u3c_.n23() - n22c - rec.q1;

      if (k1 != k2)
        return {false, false, 0};

      if ((k1 < rec.kmin) || (k1 > rec.kmax))
        return {false, false, 0};

      assert(rec.first_row + k1 < coeffs_.rows());

      std::size_t row = rec.first_row + k1 - rec.kmin;

      return {true, computed_[row] == rhomax_, row};
    }

    double get(const Accessor& acc, int rho) const
    {
      assert(acc.allowed == true);
      assert(acc.computed == true);

      return coeffs_(acc.row, rho);
    }

    void set(const Accessor& acc, int rho, double value)
    {
      assert(acc.allowed == true);

      coeffs_(acc.row, rho) = value;

      computed_[acc.row]++;
    }

    auto get_row(const Accessor& acc) const
    {
      assert(acc.allowed == true);
      assert(acc.computed == true);

      return coeffs_.row(acc.row);
    }

    template <typename T>
    void set_row(const Accessor& acc, T&& values)
    {
      assert(acc.allowed == true);

      coeffs_.row(acc.row) = std::forward<T>(values);

      computed_[acc.row] = rhomax_;
    }

    std::size_t size() const { return coeffs_.rows(); }
  };

  // expansion of a coupled gelfand basis state (u3c) in terms of those
  // uncoupled ones based on the weight selection rule; hence, this expansion is
  // the same for those other coupled states that have the same weight with the
  // state being expanded. Here indices of the uncoupled states are stored.
  std::vector<std::pair<int, int>>
  expand_to_uncoupled(const u3::gelfand& g3,
                      const std::vector<u3::gelfand>& gArr1,
                      const std::vector<u3::gelfand>& gArr2)
  {
    std::vector<std::pair<int, int>> indices;
    for (int i = 0; i < gArr1.size(); i++)
    {
      for (int j = 0; j < gArr2.size(); j++)
      {
        if ((gArr1[i].e22() + gArr2[j].e22()) == g3.e22() &&
            (gArr1[i].e11() + gArr2[j].e11()) == g3.e11())
        {
          indices.emplace_back(i, j);
        }
      }
    }
    return indices;
  }

  // build the matrix representating a system of linear equations resulting from
  // the action two generators E23 and E12 that annihiliate the highest weight
  // state
  auto E23E12(const std::vector<u3::gelfand>& gArr1,
              const std::vector<u3::gelfand>& gArr2,
              std::vector<std::pair<int, int>>& indices)
  {
    int col_num = indices.size();

    std::vector<std::pair<int, int>> new_indices; // store the parent states

    for (int k = 0; k < indices.size(); k++)
    {
      int g1 = indices[k].first, g2 = indices[k].second;

      // parent of the first state
      u3::gelfand g1_n11p1 = gArr1[g1];
      g1_n11p1.n11(g1_n11p1.n11() + 1);
      if (g1_n11p1.check())
      {
        auto i = u3::findGelfand(g1_n11p1, gArr1);
        std::pair<int, int> p(i, g2);
        auto it = std::find(new_indices.begin(), new_indices.end(), p);
        if (it == new_indices.end())
          new_indices.push_back(p);
      }

      u3::gelfand g1_n12p1 = gArr1[g1];
      g1_n12p1.n12(g1_n12p1.n12() + 1);
      if (g1_n12p1.check())
      {
        auto i = u3::findGelfand(g1_n12p1, gArr1);
        std::pair<int, int> p(i, g2);
        auto it = std::find(new_indices.begin(), new_indices.end(), p);
        if (it == new_indices.end())
          new_indices.push_back(p);
      }

      u3::gelfand g1_n22p1 = gArr1[g1];
      g1_n22p1.n22(g1_n22p1.n22() + 1);
      if (g1_n22p1.check())
      {
        auto i = u3::findGelfand(g1_n22p1, gArr1);
        std::pair<int, int> p(i, g2);
        auto it = std::find(new_indices.begin(), new_indices.end(), p);
        if (it == new_indices.end())
          new_indices.push_back(p);
      }

      // parent of the second state
      u3::gelfand g2_n11p1 = gArr2[g2];
      g2_n11p1.n11(g2_n11p1.n11() + 1);
      if (g2_n11p1.check())
      {
        auto i = u3::findGelfand(g2_n11p1, gArr2);
        std::pair<int, int> p(g1, i);
        auto it = std::find(new_indices.begin(), new_indices.end(), p);
        if (it == new_indices.end())
          new_indices.push_back(p);
      }

      u3::gelfand g2_n12p1 = gArr2[g2];
      g2_n12p1.n12(g2_n12p1.n12() + 1);
      if (g2_n12p1.check())
      {
        auto i = u3::findGelfand(g2_n12p1, gArr2);
        std::pair<int, int> p(g1, i);
        auto it = std::find(new_indices.begin(), new_indices.end(), p);
        if (it == new_indices.end())
          new_indices.push_back(p);
      }

      u3::gelfand g2_n22p1 = gArr2[g2];
      g2_n22p1.n22(g2_n22p1.n22() + 1);
      if (g2_n22p1.check())
      {
        auto i = u3::findGelfand(g2_n22p1, gArr2);
        std::pair<int, int> p(g1, i);
        auto it = std::find(new_indices.begin(), new_indices.end(), p);
        if (it == new_indices.end())
          new_indices.push_back(p);
      }
    }

    // for (const auto& p : new_indices) {
    // std::cout << "(" << p.first << ", " << p.second << ")\n";
    //}

    int row_num = new_indices.size();

#ifdef U3LIB_SPARSE_MAT
    std::vector<Eigen::Triplet<double>> ija;
#else
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(row_num, col_num);
#endif

    for (int i = 0; i < row_num; i++)
    {
      int new_g1 = new_indices[i].first, new_g2 = new_indices[i].second;
      for (int j = 0; j < col_num; j++)
      {
        int g1 = indices[j].first, g2 = indices[j].second;
        // action of E23
        if (u3::E23_12(gArr1[g1], gArr1[new_g1]) && new_g2 == g2)
#ifdef U3LIB_SPARSE_MAT
          ija.emplace_back(i, j, gArr1[g1].e23_12());
#else
          mat(i, j) = gArr1[g1].e23_12();
#endif
        else if (u3::E23_22(gArr1[g1], gArr1[new_g1]) && new_g2 == g2)
#ifdef U3LIB_SPARSE_MAT
          ija.emplace_back(i, j, gArr1[g1].e23_22());
#else
          mat(i, j) = gArr1[g1].e23_22();
#endif
        else if (new_g1 == g1 && u3::E23_12(gArr2[g2], gArr2[new_g2]))
#ifdef U3LIB_SPARSE_MAT
          ija.emplace_back(i, j, gArr2[g2].e23_12());
#else
          mat(i, j) = gArr2[g2].e23_12();
#endif
        else if (new_g1 == g1 && u3::E23_22(gArr2[g2], gArr2[new_g2]))
#ifdef U3LIB_SPARSE_MAT
          ija.emplace_back(i, j, gArr2[g2].e23_22());
#else
          mat(i, j) = gArr2[g2].e23_22();
#endif

        // action of E12
        else if (u3::E12(gArr1[g1], gArr1[new_g1]) && new_g2 == g2)
#ifdef U3LIB_SPARSE_MAT
          ija.emplace_back(i, j, gArr1[g1].e12());
#else
          mat(i, j) = gArr1[g1].e12();
#endif
        else if (new_g1 == g1 && u3::E12(gArr2[g2], gArr2[new_g2]))
#ifdef U3LIB_SPARSE_MAT
          ija.emplace_back(i, j, gArr2[g2].e12());
#else
          mat(i, j) = gArr2[g2].e12();
#endif

        else
#ifdef U3LIB_SPARSE_MAT
          ;
#else
          mat(i, j) = 0.0;
#endif
      }
    }

#ifdef U3LIB_SPARSE_MAT
    Eigen::SparseMatrix<double> mat(row_num, col_num);
    mat.setFromTriplets(ija.begin(), ija.end());
#endif

    return mat;
  }

  // we only want to store reduced U(3)->U(2) coupling coefficients,
  // so the label n11 does not matter, thus we store the reduced coefficients
  // in the coeffs matrix at row = (ia*dim(u3b) + ib)*dim(u3c) + ic
  // where ia, ib, ic are positions of U(2) highest wieght Gelfand states
  // defined as
  // | n13 n23 n33 \
  // |  n12  n22    >
  // |     n12     /
  // Given a generic Gelfand state
  // | n13 n23 n33 \
  // |  n12  n22    >
  // |     n11     /
  // at position i, the position of the U(2) highest weight state
  // can be traced back to i_hw = i - (n12 - n11) = i - n12 + n11

  // highest-weight reduced coefficients
  // / n13a n23a n33a ; n13b n23b n33b || n13c n23c n33c \
  // \   n12a  n22a       n12b  n22b   ||   n13c n23c    /
  void HighestWeightReducedCoeffs(const std::vector<u3::gelfand>& gArr1,
                                  const std::vector<u3::gelfand>& gArr2,
                                  const std::vector<u3::gelfand>& gArr3,
                                  u3::CGC_SO3& cgcs, WigU3U2_coeffs& coeffs)
  {
    std::vector<std::pair<int, int>> indices =
        expand_to_uncoupled(gArr3[0], gArr1, gArr2);

    auto mat = E23E12(gArr1, gArr2, indices);

#ifdef U3LIB_TIME_LOG
    auto t0 = std::chrono::high_resolution_clock::now();
#endif

#ifdef U3LIB_NULLSPACE_LOG
    std::cout << "mat: " << mat.rows() << " x " << mat.cols() << std::endl;
    std::cout << "\telements: " << mat.rows() * mat.cols() << std::endl;

#ifdef U3LIB_SPARSE_MAT
    long nnz = mat.nonZeros();
#else
    long nnz = 0;
    for (auto x : mat.reshaped())
      if (x != 0.0)
        nnz++;
#endif
    std::cout << "\tnonzero elements: " << nnz << std::endl;
#endif

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nullspace =
        null_space(mat);

#ifdef U3LIB_NULLSPACE_LOG
    std::cout << "nullspace: " << nullspace.rows() << " x " << nullspace.cols() << std::endl;
    std::cout << "\telements: " << nullspace.rows() * nullspace.cols() << std::endl;
    nnz = 0;
    for (auto x : nullspace.reshaped())
      if (x != 0.0)
        nnz++;
    std::cout << "\tnonzero elements: " << nnz << std::endl;
#endif

#ifdef U3LIB_TIME_LOG
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = t1 - t0;
    std::cerr << "null_space time: " << duration.count() << " [s]" << std::endl;
#endif

    for (int i = 0; i < indices.size(); i++)
    {
      int i1 = indices[i].first - gArr1[indices[i].first].n12() +
               gArr1[indices[i].first].n11();
      int i2 = indices[i].second - gArr2[indices[i].second].n12() +
               gArr2[indices[i].second].n11();

      const auto& gA = gArr1[i1];
      const auto& gB = gArr2[i2];
      const auto& gC = gArr3[0];

      auto acc = coeffs.access(gA.n12(), gA.n22(), gB.n12(), gB.n22(), gC.n12(),
                               gC.n22());

      if (acc.computed == false)
      {
        int jj1 = gArr1[indices[i].first].Lambdax2(),
            jj2 = gArr2[indices[i].second].Lambdax2(),
            jj3 = gArr3[0].Lambdax2(),
            mm1 = gArr1[indices[i].first].MLambdax2(),
            mm2 = gArr2[indices[i].second].MLambdax2(),
            mm3 = gArr3[0].MLambdax2();

        if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
            (jj1 + jj2) >= jj3)
        {
          double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

          assert(cgc != 0.0);

          coeffs.set_row(acc, nullspace.row(i) / cgc);
        }
      }
    }

    return;
  }

  // lower-weight reduced coefficients
  // / n13a n23a n33a ; n13b n23b n33b || n13c n23c n33c \
  // \   n12a  n22a       n12b  n22b   ||   n12c n22c    /

  // build with E32
  void E32(int parent_g3, int& row, int rhomax, 
           std::vector<int>& same_weight_gArr3, 
           std::vector<std::pair<int, int>>& uncoupled_same_weight_gArr3,
           const std::vector<u3::gelfand>& gArr1,
           const std::vector<u3::gelfand>& gArr2,
           const std::vector<u3::gelfand>& gArr3,
           Eigen::MatrixXd& b, u3::CGC_SO3& cgcs, 
           std::vector<Eigen::MatrixXd>& E32E31, 
           WigU3U2_coeffs& coeffs)
  {
    int jj3 = gArr3[parent_g3].Lambdax2(),
        mm3 = gArr3[parent_g3].MLambdax2();

    // loop over daughter states
    for (int j = 0; j < same_weight_gArr3.size(); j++)
      if (u3::E23_12(gArr3[same_weight_gArr3[j]], gArr3[parent_g3]))
        b(row, j) = gArr3[parent_g3].e32_12();
      else if (u3::E23_22(gArr3[same_weight_gArr3[j]], gArr3[parent_g3]))
        b(row, j) = gArr3[parent_g3].e32_22();
      else
        b(row, j) = 0.0;

    for (int j = 0; j < uncoupled_same_weight_gArr3.size(); j++)
    {
      int g1 = uncoupled_same_weight_gArr3[j].first,
          g2 = uncoupled_same_weight_gArr3[j].second;

      // act on the first state
      u3::gelfand g1_n12p1 = gArr1[g1];
      g1_n12p1.n12(g1_n12p1.n12() + 1);
      if (g1_n12p1.check())
      {
        int jj1 = g1_n12p1.Lambdax2(), mm1 = g1_n12p1.MLambdax2(),
            jj2 = gArr2[g2].Lambdax2(), mm2 = gArr2[g2].MLambdax2();

        if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
            (jj1 + jj2) >= jj3)
        {
          double factor = g1_n12p1.e32_12();

          double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

          for (int r = 0; r < rhomax; r++)
          {
            auto acc = coeffs.access(g1_n12p1.n12(), g1_n12p1.n22(),
                                     gArr2[g2].n12(), gArr2[g2].n22(),
                                     gArr3[parent_g3].n12(),
                                     gArr3[parent_g3].n22());

            E32E31[r](row, j) += coeffs.get(acc, r) * cgc * factor;
          }
        }
      }

      u3::gelfand g1_n22p1 = gArr1[g1];
      g1_n22p1.n22(g1_n22p1.n22() + 1);
      if (g1_n22p1.check())
      {
        int jj1 = g1_n22p1.Lambdax2(), mm1 = g1_n22p1.MLambdax2(),
            jj2 = gArr2[g2].Lambdax2(), mm2 = gArr2[g2].MLambdax2();

        if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
            (jj1 + jj2) >= jj3)
        {
          double factor = g1_n22p1.e32_22();

          double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

          for (int r = 0; r < rhomax; r++)
          {
            auto acc = coeffs.access(g1_n22p1.n12(), g1_n22p1.n22(),
                                     gArr2[g2].n12(), gArr2[g2].n22(),
                                     gArr3[parent_g3].n12(),
                                     gArr3[parent_g3].n22());

            E32E31[r](row, j) += coeffs.get(acc, r) * cgc * factor;
          }
        }
      }

      // act on the second state
      u3::gelfand g2_n12p1 = gArr2[g2];
      g2_n12p1.n12(g2_n12p1.n12() + 1);
      if (g2_n12p1.check())
      {
        int jj2 = g2_n12p1.Lambdax2(), mm2 = g2_n12p1.MLambdax2(),
            jj1 = gArr1[g1].Lambdax2(), mm1 = gArr1[g1].MLambdax2();

        if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
            (jj1 + jj2) >= jj3)
        {
          double factor = g2_n12p1.e32_12();

          double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

          for (int r = 0; r < rhomax; r++)
          {
            auto acc = coeffs.access(gArr1[g1].n12(), gArr1[g1].n22(),
                                     g2_n12p1.n12(), g2_n12p1.n22(),
                                     gArr3[parent_g3].n12(),
                                     gArr3[parent_g3].n22());

            E32E31[r](row, j) += coeffs.get(acc, r) * cgc * factor;
          }
        }
      }

      u3::gelfand g2_n22p1 = gArr2[g2];
      g2_n22p1.n22(g2_n22p1.n22() + 1);
      if (g2_n22p1.check())
      {
        int jj2 = g2_n22p1.Lambdax2(), mm2 = g2_n22p1.MLambdax2(),
            jj1 = gArr1[g1].Lambdax2(), mm1 = gArr1[g1].MLambdax2();

        if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
            (jj1 + jj2) >= jj3)
        {
          double factor = g2_n22p1.e32_22();

          double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

          for (int r = 0; r < rhomax; r++)
          {
            auto acc = coeffs.access(gArr1[g1].n12(), gArr1[g1].n22(),
                                     g2_n22p1.n12(), g2_n22p1.n22(),
                                     gArr3[parent_g3].n12(),
                                     gArr3[parent_g3].n22());

            E32E31[r](row, j) += coeffs.get(acc, r) * cgc * factor;
          }
        }
      }
    }
    row++;
    return;
  }

  // build with E31
  void E31(int parent_g3, int& row, int rhomax, 
           std::vector<int>& same_weight_gArr3, 
           std::vector<std::pair<int, int>>& uncoupled_same_weight_gArr3,
           const std::vector<u3::gelfand>& gArr1,
           const std::vector<u3::gelfand>& gArr2,
           const std::vector<u3::gelfand>& gArr3,
           Eigen::MatrixXd& b, u3::CGC_SO3& cgcs, 
           std::vector<Eigen::MatrixXd>& E32E31, 
           WigU3U2_coeffs& coeffs)
  {
    int jj3 = gArr3[parent_g3].Lambdax2(),
        mm3 = gArr3[parent_g3].MLambdax2();

    // loop over daughter states
    for (int j = 0; j < same_weight_gArr3.size(); j++)
      if (u3::E13_12(gArr3[same_weight_gArr3[j]], gArr3[parent_g3]))
        b(row, j) = gArr3[parent_g3].e31_12();
      else if (u3::E13_22(gArr3[same_weight_gArr3[j]], gArr3[parent_g3]))
        b(row, j) = gArr3[parent_g3].e31_22();
      else
        b(row, j) = 0.0;

    for (int j = 0; j < uncoupled_same_weight_gArr3.size(); j++)
    {
      int g1 = uncoupled_same_weight_gArr3[j].first,
          g2 = uncoupled_same_weight_gArr3[j].second;

      // act on the first state
      u3::gelfand g1_n12n11p1 = gArr1[g1];
      g1_n12n11p1.n12(g1_n12n11p1.n12() + 1);
      g1_n12n11p1.n11(g1_n12n11p1.n11() + 1);
      if (g1_n12n11p1.check())
      {
        int jj1 = g1_n12n11p1.Lambdax2(), mm1 = g1_n12n11p1.MLambdax2(),
            jj2 = gArr2[g2].Lambdax2(), mm2 = gArr2[g2].MLambdax2();

        if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
            (jj1 + jj2) >= jj3)
        {
          double factor = g1_n12n11p1.e31_12();

          double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

          for (int r = 0; r < rhomax; r++)
          {
            auto acc = coeffs.access(g1_n12n11p1.n12(), g1_n12n11p1.n22(),
                                     gArr2[g2].n12(), gArr2[g2].n22(),
                                     gArr3[parent_g3].n12(),
                                     gArr3[parent_g3].n22());

            E32E31[r](row, j) += coeffs.get(acc, r) * cgc * factor;
          }
        }
      }

      u3::gelfand g1_n22n11p1 = gArr1[g1];
      g1_n22n11p1.n22(g1_n22n11p1.n22() + 1);
      g1_n22n11p1.n11(g1_n22n11p1.n11() + 1);
      if (g1_n22n11p1.check())
      {
        int jj1 = g1_n22n11p1.Lambdax2(), mm1 = g1_n22n11p1.MLambdax2(),
            jj2 = gArr2[g2].Lambdax2(), mm2 = gArr2[g2].MLambdax2();

        if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
            (jj1 + jj2) >= jj3)
        {
          double factor = g1_n22n11p1.e31_22();

          double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

          for (int r = 0; r < rhomax; r++)
          {
            auto acc = coeffs.access(g1_n22n11p1.n12(), g1_n22n11p1.n22(),
                                     gArr2[g2].n12(), gArr2[g2].n22(),
                                     gArr3[parent_g3].n12(),
                                     gArr3[parent_g3].n22());

            E32E31[r](row, j) += coeffs.get(acc, r) * cgc * factor;
          }
        }
      }

      // act on the second state
      u3::gelfand g2_n12n11p1 = gArr2[g2];
      g2_n12n11p1.n12(g2_n12n11p1.n12() + 1);
      g2_n12n11p1.n11(g2_n12n11p1.n11() + 1);
      if (g2_n12n11p1.check())
      {
        int jj2 = g2_n12n11p1.Lambdax2(), mm2 = g2_n12n11p1.MLambdax2(),
            jj1 = gArr1[g1].Lambdax2(), mm1 = gArr1[g1].MLambdax2();

        if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
            (jj1 + jj2) >= jj3)
        {
          double factor = g2_n12n11p1.e31_12();

          double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

          for (int r = 0; r < rhomax; r++)
          {
            auto acc = coeffs.access(gArr1[g1].n12(), gArr1[g1].n22(),
                                     g2_n12n11p1.n12(), g2_n12n11p1.n22(),
                                     gArr3[parent_g3].n12(),
                                     gArr3[parent_g3].n22());

            E32E31[r](row, j) += coeffs.get(acc, r) * cgc * factor;
          }
        }
      }

      u3::gelfand g2_n22n11p1 = gArr2[g2];
      g2_n22n11p1.n22(g2_n22n11p1.n22() + 1);
      g2_n22n11p1.n11(g2_n22n11p1.n11() + 1);
      if (g2_n22n11p1.check())
      {
        int jj2 = g2_n22n11p1.Lambdax2(), mm2 = g2_n22n11p1.MLambdax2(),
            jj1 = gArr1[g1].Lambdax2(), mm1 = gArr1[g1].MLambdax2();

        if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
            (jj1 + jj2) >= jj3)
        {
          double factor = g2_n22n11p1.e31_22();

          double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

          for (int r = 0; r < rhomax; r++)
          {
            auto acc = coeffs.access(gArr1[g1].n12(), gArr1[g1].n22(),
                                     g2_n22n11p1.n12(), g2_n22n11p1.n22(),
                                     gArr3[parent_g3].n12(),
                                     gArr3[parent_g3].n22());

            E32E31[r](row, j) += coeffs.get(acc, r) * cgc * factor;
          }
        }
      }
    }
    row++;
    return;
  }

  void LowerWeightReducedCoeffs(const std::vector<u3::gelfand>& gArr1,
                                const std::vector<u3::gelfand>& gArr2,
                                const std::vector<u3::gelfand>& gArr3,
                                int rhomax, std::vector<int>& same_weight_gArr3,
                                u3::CGC_SO3& cgcs, WigU3U2_coeffs& coeffs)
  {
    for (int i = 0; i < same_weight_gArr3.size(); i++)
    {
      same_weight_gArr3[i] += gArr3[same_weight_gArr3[i]].n11() -
                              gArr3[same_weight_gArr3.back()].n11();
    }

    if (!same_weight_gArr3.empty())
    {
      // the coupled states having the same weight have the same expansion of
      // uncoupled states, whose indices are recorded
      std::vector<std::pair<int, int>> uncoupled_same_weight_gArr3 =
          expand_to_uncoupled(gArr3[same_weight_gArr3[0]], gArr1, gArr2);

      // matrix b is independent of outer multiplicity
      Eigen::MatrixXd b(same_weight_gArr3.size(), same_weight_gArr3.size());

      // the lowering matrix from the action of E32 and E31, which depends on
      // outer multiplicity
      std::vector<Eigen::MatrixXd> E32E31(rhomax);

      // set matrix elements to zero
      for (int r = 0; r < rhomax; r++)
      {
        E32E31[r] = Eigen::MatrixXd(same_weight_gArr3.size(),
                                    uncoupled_same_weight_gArr3.size());
        E32E31[r].setZero();
      }

      // find parent states and operators to get to lower weight states
      // (daughter states)
      std::vector<std::pair<int, int>> actOnParent;

      int row = 0;

      // loop over daughter states
      for (int i = 0; i < same_weight_gArr3.size(); i++)
      {
        u3::gelfand g3_n12p1 = gArr3[same_weight_gArr3[i]];
        g3_n12p1.n12(g3_n12p1.n12() + 1);
        if (g3_n12p1.check())
        {
          auto j = u3::findGelfand(g3_n12p1, gArr3);
          std::pair<int, int> p(0, j);
          auto it = std::find(actOnParent.begin(), actOnParent.end(), p);
          if (it == actOnParent.end())
          {
            actOnParent.push_back(p);
            E32(j, row, rhomax, same_weight_gArr3, uncoupled_same_weight_gArr3,
                gArr1, gArr2, gArr3, b, cgcs, E32E31, coeffs);
            continue;
          }
        }

        u3::gelfand g3_n22p1 = gArr3[same_weight_gArr3[i]];
        g3_n22p1.n22(g3_n22p1.n22() + 1);
        if (g3_n22p1.check())
        {
          auto j = u3::findGelfand(g3_n22p1, gArr3);
          std::pair<int, int> p(0, j);
          auto it = std::find(actOnParent.begin(), actOnParent.end(), p);
          if (it == actOnParent.end())
          {
            actOnParent.push_back(p);
            E32(j, row, rhomax, same_weight_gArr3, uncoupled_same_weight_gArr3,
                gArr1, gArr2, gArr3, b, cgcs, E32E31, coeffs);
            continue;
          }
        }

        u3::gelfand g3_n12n11p1 = gArr3[same_weight_gArr3[i]];
        g3_n12n11p1.n12(g3_n12n11p1.n12() + 1);
        g3_n12n11p1.n11(g3_n12n11p1.n11() + 1);
        if (g3_n12n11p1.check())
        {
          auto j = u3::findGelfand(g3_n12n11p1, gArr3);
          std::pair<int, int> p(1, j);
          auto it = std::find(actOnParent.begin(), actOnParent.end(), p);
          if (it == actOnParent.end())
          {
            actOnParent.push_back(p);
            E31(j, row, rhomax, same_weight_gArr3, uncoupled_same_weight_gArr3,
                gArr1, gArr2, gArr3, b, cgcs, E32E31, coeffs);
            continue;
          }
        }

        u3::gelfand g3_n22n11p1 = gArr3[same_weight_gArr3[i]];
        g3_n22n11p1.n22(g3_n22n11p1.n22() + 1);
        g3_n22n11p1.n11(g3_n22n11p1.n11() + 1);
        if (g3_n22n11p1.check())
        {
          auto j = u3::findGelfand(g3_n22n11p1, gArr3);
          std::pair<int, int> p(1, j);
          auto it = std::find(actOnParent.begin(), actOnParent.end(), p);
          if (it == actOnParent.end())
          {
            actOnParent.push_back(p);
            E31(j, row, rhomax, same_weight_gArr3, uncoupled_same_weight_gArr3,
                gArr1, gArr2, gArr3, b, cgcs, E32E31, coeffs);
            continue;
          }
        }
      }

      // std::cout << b.rows() << "," << b.cols() << std::endl;

      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> factorize;
      factorize.compute(b);

      // Eigen::MatrixXd b_inverse = inverse(b);

      std::vector<Eigen::MatrixXd> coeffs_temp(rhomax);

      for (int r = 0; r < rhomax; r++)
        coeffs_temp[r] = factorize.solve(E32E31[r]);
      // coeffs_temp[r] = b_inverse * E32E31[r];

      Eigen::RowVectorXd coeffs_row(rhomax); // define here to avoid allocation

      for (int i = 0; i < same_weight_gArr3.size(); i++)
      {
        for (int j = 0; j < uncoupled_same_weight_gArr3.size(); j++)
        {

          int i1 = uncoupled_same_weight_gArr3[j].first,
              i2 = uncoupled_same_weight_gArr3[j].second,
              i3 = same_weight_gArr3[i];

          int i1_hw = i1 - gArr1[i1].n12() + gArr1[i1].n11(),
              i2_hw = i2 - gArr2[i2].n12() + gArr2[i2].n11(),
              i3_hw = i3 - gArr3[i3].n12() + gArr3[i3].n11();

          const auto& gA = gArr1[i1_hw];
          const auto& gB = gArr2[i2_hw];
          const auto& gC = gArr3[i3_hw];

          auto acc = coeffs.access(gA.n12(), gA.n22(), gB.n12(), gB.n22(),
                                   gC.n12(), gC.n22());

          if (acc.computed == false)
          {
            int jj1 = gArr1[i1].Lambdax2(), jj2 = gArr2[i2].Lambdax2(),
                jj3 = gArr3[i3].Lambdax2(), mm1 = gArr1[i1].MLambdax2(),
                mm2 = gArr2[i2].MLambdax2(), mm3 = gArr3[i3].MLambdax2();

            if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
                (jj1 + jj2) >= jj3)
            {
              double cgc = cgcs.compute(jj1, jj2, jj3, mm1, mm2, mm3);

              for (int r = 0; r < rhomax; r++)
              {
                if (std::fabs(coeffs_temp[r](i, j)) < 1e-10)
                  coeffs_temp[r](i, j) = 0.0;

                assert(cgc != 0.0);

                coeffs_row(r) = coeffs_temp[r](i, j) / cgc;
              }

              coeffs.set_row(acc, coeffs_row);
            }
          }
        }
      }
    }

    return;
  }

  // all coupling coefficients
  auto WigU3U2(const u3::irrep& u3A, const u3::irrep& u3B, u3::irrep& u3C,
               int rhomax, int w2min, std::vector<u3::gelfand>& gArrA,
               std::vector<u3::gelfand>& gArrB, std::vector<u3::gelfand>& gArrC)
  {
    WigU3U2_coeffs coeffs(u3A, u3B, u3C, rhomax);

    int w2max = u3C.n13() + u3C.n23();
    if (w2min == -1) // calculate all
      w2min = u3C.n23() + u3C.n33();

    std::vector<std::vector<int>> all_same_weight_gArrC(w2max - w2min + 1);
    for (int i = 0; i < gArrC.size(); i++)
      if (gArrC[i].n12() == gArrC[i].n11())
        all_same_weight_gArrC[w2max - gArrC[i].n12() - gArrC[i].n22()]
            .push_back(i);

    u3::CGC_SO3 cgcs;

    u3::HighestWeightReducedCoeffs(gArrA, gArrB, gArrC, cgcs, coeffs);

#ifdef U3LIB_TIME_LOG
    auto t0 = std::chrono::high_resolution_clock::now();
#endif

    for (int w2 = w2max - 1; w2 >= w2min; w2--)
    {
      u3::LowerWeightReducedCoeffs(gArrA, gArrB, gArrC, rhomax,
                                   all_same_weight_gArrC[w2max - w2], cgcs,
                                   coeffs);
    }

#ifdef U3LIB_TIME_LOG
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = t1 - t0;
    std::cerr << "LowerWeightReducedCoeffs total time: " << duration.count()
              << " [s]" << std::endl;
#endif

    u3C.normalize();
    for (int i = 0; i < gArrC.size(); i++)
      gArrC[i].normalize();

    // once bases are normalized, coeffs needs to be as well:
    coeffs.normalize();

    return coeffs;
  }

} // namespace u3

#endif

/*
    // Initialize max values
    int maxFirst = indices.back().first;
    int maxSecond = std::numeric_limits<int>::min();

    // Iterate through the vector to find the maximum
    for (const auto& p : indices)
      maxSecond = std::max(maxSecond, p.second); // Update maxSecond

    for (int i = 0; i < maxFirst; i++)
    {
      for (int j = 0; j < maxSecond; j++)
      {
        for (int k = 0; k < indices.size(); k++)
        {
          int g1 = indices[k].first, g2 = indices[k].second;

          // E23 and E12 act on the first state
          if ((u3::E23_12(gArr1[g1], gArr1[i]) && g2 == j) ||
              (u3::E23_22(gArr1[g1], gArr1[i]) && g2 == j) ||
              (u3::E12(gArr1[g1], gArr1[i]) && g2 == j))
            new_indices.emplace_back(i, g2);

          // E23 and E12 act on the second state
          else if ((u3::E23_12(gArr2[g2], gArr2[j]) && g1 == i) ||
                   (u3::E23_22(gArr2[g2], gArr2[j]) && g1 == i) ||
                   (u3::E12(gArr2[g2], gArr2[j]) && g1 == i))
            new_indices.emplace_back(g1, j);
        }
      }
    }
*/

/*
      // loop over states that may be parent
      for (int j = 0; j < same_weight_gArr3.back(); j++)
      {
        // loop over daughter states
        for (int i = 0; i < same_weight_gArr3.size(); i++)
        {
          if (u3::E23_12(gArr3[same_weight_gArr3[i]], gArr3[j]) ||
              u3::E23_22(gArr3[same_weight_gArr3[i]], gArr3[j]))
          {
            std::pair p(0, j);
            auto it = std::find(actOnParent.begin(), actOnParent.end(), p);
            if (it == actOnParent.end())
              actOnParent.push_back(p);
          }

          else if (u3::E13_12(gArr3[same_weight_gArr3[i]], gArr3[j]) ||
                   u3::E13_22(gArr3[same_weight_gArr3[i]], gArr3[j]))
          {
            std::pair p(1, j);
            auto it = std::find(actOnParent.begin(), actOnParent.end(), p);
            if (it == actOnParent.end())
              actOnParent.push_back(p);
          }
        }
      }
*/
