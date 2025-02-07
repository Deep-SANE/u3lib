#ifndef U3LIB_U3U2U1XU3U2U1_H
#define U3LIB_U3U2U1XU3U2U1_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "decomposition.h"
#include "U3basis.h"
#include "utils/inverseSolver.h"
#include "utils/nullspaceSolver.h"

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
          std::pair<int, int> p(i, j);
          indices.push_back(p);
          // std::cout << "-" << i << "," << j << std::endl;
        }
      }
    }
    return indices;
  }

  // build the matrix representating a system of linear equations resulting from
  // the action two generators E23 and E12 that annihiliate the highest weight
  // state
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  E23E12(const std::vector<u3::gelfand>& gArr1,
         const std::vector<u3::gelfand>& gArr2,
         std::vector<std::pair<int, int>>& indices)
  {
    int col_num = indices.size();

    std::vector<std::pair<int, int>> new_indices; // store the parent states

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

    int row_num = new_indices.size();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat(row_num, col_num);

    for (int i = 0; i < row_num; i++)
    {
      // std::cout << new_indices[i].first << "," << new_indices[i].second <<
      // std::endl;

      int new_g1 = new_indices[i].first, new_g2 = new_indices[i].second;
      for (int j = 0; j < col_num; j++)
      {
        int g1 = indices[j].first, g2 = indices[j].second;
        // action of E23
        if (u3::E23_12(gArr1[g1], gArr1[new_g1]) && new_g2 == g2)
          mat(i, j) = gArr1[g1].e23_12();
        else if (u3::E23_22(gArr1[g1], gArr1[new_g1]) && new_g2 == g2)
          mat(i, j) = gArr1[g1].e23_22();
        else if (new_g1 == g1 && u3::E23_12(gArr2[g2], gArr2[new_g2]))
          mat(i, j) = gArr2[g2].e23_12();
        else if (new_g1 == g1 && u3::E23_22(gArr2[g2], gArr2[new_g2]))
          mat(i, j) = gArr2[g2].e23_22();

        // action of E12
        else if (u3::E12(gArr1[g1], gArr1[new_g1]) && new_g2 == g2)
          mat(i, j) = gArr1[g1].e12();
        else if (new_g1 == g1 && u3::E12(gArr2[g2], gArr2[new_g2]))
          mat(i, j) = gArr2[g2].e12();

        else
          mat(i, j) = 0.0;
      }
    }

    return mat;
  }

  // highest-weight coefficients < G1; G2| HW3 >
  void HighestWeightCoeffs(
      const std::vector<u3::gelfand>& gArr1,
      const std::vector<u3::gelfand>& gArr2,
      const std::vector<u3::gelfand>& gArr3,
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& coeffs)
  {

    std::vector<std::pair<int, int>> indices =
        expand_to_uncoupled(gArr3[0], gArr1, gArr2);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat =
        E23E12(gArr1, gArr2, indices);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nullspace =
        null_space(mat);

    for (int i = 0; i < indices.size(); i++)
      coeffs.row((indices[i].first * gArr2.size() + indices[i].second) *
                 gArr3.size()) = nullspace.row(i);

    return;
  }

  // lower-weight coefficients < G1; G2 | G3 >

  void LowerWeightCoeffs(
      int w3, int w1, const std::vector<u3::gelfand>& gArr1,
      const std::vector<u3::gelfand>& gArr2,
      const std::vector<u3::gelfand>& gArr3,
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& coeffs)
  {

    /// record indices of all coupled states that have the same weight
    std::vector<int> same_weight_gArr3;

    for (int i = 0; i < gArr3.size(); i++)
      if (gArr3[i].e11() == w1 && gArr3[i].e33() == w3)
        same_weight_gArr3.push_back(i);

    if (same_weight_gArr3.size() != 0)
    {
      // the coupled states having the same weight have the same expansion of
      // uncoupled states, whose indices are recorded
      std::vector<std::pair<int, int>> uncoupled_same_weight_gArr3 =
          expand_to_uncoupled(gArr3[same_weight_gArr3[0]], gArr1, gArr2);

      // matrix b is independent of outer multiplicity
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b(
          same_weight_gArr3.size(), same_weight_gArr3.size());

      // the lowering matrix from the action of E21 and E32, which depends on
      // outer multiplicity
      int rmax = coeffs.cols();

      std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> E21E32(
          rmax);

      // set matrix elements to zero
      for (int r = 0; r < rmax; r++)
      {
        E21E32[r] = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(
            same_weight_gArr3.size(), uncoupled_same_weight_gArr3.size());
        E21E32[r].setZero();
      }

      for (int i = 0; i < same_weight_gArr3.size(); i++)
      {
        int flag = 0;

        // build with generator E21

        u3::gelfand g3_n11p1 = gArr3[same_weight_gArr3[i]];
        g3_n11p1.n11(g3_n11p1.n11() + 1); // this is a parent state

        if (g3_n11p1.check())
        {
          std::vector<std::pair<int, int>> uncoupled_g3n11p1 =
              expand_to_uncoupled(g3_n11p1, gArr1, gArr2);

          for (int j = 0; j < same_weight_gArr3.size(); j++)
            if (u3::E12(gArr3[same_weight_gArr3[j]], g3_n11p1))
              b(i, j) = g3_n11p1.e21();
            else
              b(i, j) = 0.0;

          int k = u3::findGelfand(g3_n11p1, gArr3);

          for (int j = 0; j < uncoupled_same_weight_gArr3.size(); j++)
          {
            int g1 = uncoupled_same_weight_gArr3[j].first,
                g2 = uncoupled_same_weight_gArr3[j].second;

            for (int l = 0; l < uncoupled_g3n11p1.size(); l++)
            {
              int parent_g1 = uncoupled_g3n11p1[l].first,
                  parent_g2 = uncoupled_g3n11p1[l].second;

              bool cond1 =
                  (u3::E12(gArr1[g1], gArr1[parent_g1]) && g2 == parent_g2);
              bool cond2 =
                  (u3::E12(gArr2[g2], gArr2[parent_g2]) && g1 == parent_g1);

              if (cond1 || cond2)
              {
                double factor =
                    cond1 ? gArr1[parent_g1].e21() : gArr2[parent_g2].e21();

                for (int r = 0; r < rmax; r++)
                  E21E32[r](i, j) +=
                      coeffs((parent_g1 * gArr2.size() + parent_g2) *
                                     gArr3.size() +
                                 k,
                             r) *
                      factor;

                flag++;
              }
            }
          }
        }

        if (flag > 0)
          continue;

        // build with generator E32

        u3::gelfand g3_n12p1 = gArr3[same_weight_gArr3[i]];
        g3_n12p1.n12(g3_n12p1.n12() + 1); // this is a parent state

        if (g3_n12p1.check())
        {
          std::vector<std::pair<int, int>> uncoupled_g3n12p1 =
              expand_to_uncoupled(g3_n12p1, gArr1, gArr2);

          for (int j = 0; j < same_weight_gArr3.size(); j++)
            if (u3::E23_12(gArr3[same_weight_gArr3[j]], g3_n12p1))
              b(i, j) = g3_n12p1.e32_12();
            else if (u3::E23_22(gArr3[same_weight_gArr3[j]], g3_n12p1))
              b(i, j) = g3_n12p1.e32_22();
            else
              b(i, j) = 0.0;

          int k = u3::findGelfand(g3_n12p1, gArr3);

          for (int j = 0; j < uncoupled_same_weight_gArr3.size(); j++)
          {
            int g1 = uncoupled_same_weight_gArr3[j].first,
                g2 = uncoupled_same_weight_gArr3[j].second;

            for (int l = 0; l < uncoupled_g3n12p1.size(); l++)
            {
              int parent_g1 = uncoupled_g3n12p1[l].first,
                  parent_g2 = uncoupled_g3n12p1[l].second;

              bool cond1 =
                  (u3::E23_12(gArr1[g1], gArr1[parent_g1]) && g2 == parent_g2);
              bool cond2 =
                  (u3::E23_22(gArr1[g1], gArr1[parent_g1]) && g2 == parent_g2);
              bool cond3 =
                  (u3::E23_12(gArr2[g2], gArr2[parent_g2]) && g1 == parent_g1);
              bool cond4 =
                  (u3::E23_22(gArr2[g2], gArr2[parent_g2]) && g1 == parent_g1);

              if (cond1 || cond2 || cond3 || cond4)
              {
                double factor;
                if (cond1)
                  factor = gArr1[parent_g1].e32_12();
                else if (cond2)
                  factor = gArr1[parent_g1].e32_22();
                else if (cond3)
                  factor = gArr2[parent_g2].e32_12();
                else
                  gArr2[parent_g2].e32_22();

                for (int r = 0; r < rmax; r++)
                  E21E32[r](i, j) +=
                      coeffs((parent_g1 * gArr2.size() + parent_g2) *
                                     gArr3.size() +
                                 k,
                             r) *
                      factor;

                flag++;
              }
            }
          }
        }

        if (flag > 0)
          continue;

        u3::gelfand g3_n22p1 = gArr3[same_weight_gArr3[i]];
        g3_n22p1.n22(g3_n22p1.n22() + 1); // this is a parent state

        if (g3_n22p1.check())
        {
          std::vector<std::pair<int, int>> uncoupled_g3n22p1 =
              expand_to_uncoupled(g3_n22p1, gArr1, gArr2);

          for (int j = 0; j < same_weight_gArr3.size(); j++)
            if (u3::E23_12(gArr3[same_weight_gArr3[j]], g3_n22p1))
              b(i, j) = g3_n22p1.e32_12();
            else if (u3::E23_22(gArr3[same_weight_gArr3[j]], g3_n22p1))
              b(i, j) = g3_n22p1.e32_22();
            else
              b(i, j) = 0.0;

          int k = u3::findGelfand(g3_n22p1, gArr3);

          for (int j = 0; j < uncoupled_same_weight_gArr3.size(); j++)
          {
            int g1 = uncoupled_same_weight_gArr3[j].first,
                g2 = uncoupled_same_weight_gArr3[j].second;

            for (int l = 0; l < uncoupled_g3n22p1.size(); l++)
            {
              int parent_g1 = uncoupled_g3n22p1[l].first,
                  parent_g2 = uncoupled_g3n22p1[l].second;

              bool cond1 =
                  (u3::E23_12(gArr1[g1], gArr1[parent_g1]) && g2 == parent_g2);
              bool cond2 =
                  (u3::E23_22(gArr1[g1], gArr1[parent_g1]) && g2 == parent_g2);
              bool cond3 =
                  (u3::E23_12(gArr2[g2], gArr2[parent_g2]) && g1 == parent_g1);
              bool cond4 =
                  (u3::E23_22(gArr2[g2], gArr2[parent_g2]) && g1 == parent_g1);

              if (cond1 || cond2 || cond3 || cond4)
              {
                double factor;
                if (cond1)
                  factor = gArr1[parent_g1].e32_12();
                else if (cond2)
                  factor = gArr1[parent_g1].e32_22();
                else if (cond3)
                  factor = gArr2[parent_g2].e32_12();
                else
                  factor = gArr2[parent_g2].e32_22();

                for (int r = 0; r < rmax; r++)
                  E21E32[r](i, j) +=
                      coeffs((parent_g1 * gArr2.size() + parent_g2) *
                                     gArr3.size() +
                                 k,
                             r) *
                      factor;

                flag++;
              }
            }
          }
        }
      }

      // inverse b
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_inverse =
          inverse(b);

      // std::cout << "b = " << b << std::endl;

      for (int r = 0; r < rmax; r++)
      {
        // std::cout << "E21E32 = " << E21E32[r] << std::endl;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffs_temp =
            b_inverse * E21E32[r];

        for (int i = 0; i < same_weight_gArr3.size(); i++)
        {
          for (int j = 0; j < uncoupled_same_weight_gArr3.size(); j++)
          {
            if (std::fabs(coeffs_temp(i, j)) < 1e-10)
              coeffs_temp(i, j) = 0.0;

            coeffs((uncoupled_same_weight_gArr3[j].first * gArr2.size() +
                    uncoupled_same_weight_gArr3[j].second) *
                           gArr3.size() +
                       same_weight_gArr3[i],
                   r) = coeffs_temp(i, j);
          }
        }
      }
    }
  }

  // all coupling coefficients
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  WigU3U2U1(const u3::irrep& u3A, const u3::irrep& u3B, u3::irrep& u3C,
            const int& rhomax, std::vector<u3::gelfand>& gArrA,
            std::vector<u3::gelfand>& gArrB, std::vector<u3::gelfand>& gArrC)
  {
    gArrA = u3::gelfandBasis(u3A);
    gArrB = u3::gelfandBasis(u3B);
    gArrC = u3::gelfandBasis(u3C);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffs(
        gArrA.size() * gArrB.size() * gArrC.size(), rhomax);
    coeffs.setZero();

    for (int w3 = u3C.n33(); w3 <= u3C.n13(); w3++)
    {
      for (int w1 = u3C.n13(); w1 >= u3C.n33(); w1--)
        if (w3 == u3C.n33() && w1 == u3C.n13())
          u3::HighestWeightCoeffs(gArrA, gArrB, gArrC, coeffs);
        else
          u3::LowerWeightCoeffs(w3, w1, gArrA, gArrB, gArrC, coeffs);
    }

    u3C.normalize();
    for (int i = 0; i < gArrC.size(); i++)
      gArrC[i].normalize();

    return coeffs;
  }

} // namespace u3

#endif
