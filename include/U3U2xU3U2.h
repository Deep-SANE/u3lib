#ifndef U3LIB_U3U2XU3U2_H
#define U3LIB_U3U2XU3U2_H

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <tuple>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "decomposition.h"
#include "U3basis.h"
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
  // usual SO(3) Clebsh-Gordan coefficients <L1 M1; L2 M2 | L3 M3>
  double CGC_SO3(int ll1, int ll2, int ll3, int mm1, int mm2, int mm3)
  {
    if ((ll3 == (ll1 + ll2) && ((mm1 + mm2) == ll3 || (mm1 + mm2) == -ll3)) ||
        ll1 == 0 || ll2 == 0)
      return 1.0;
    else
    {
      int phase = ((-ll1 + ll2 - mm3) % 4 == 0) ? 1 : -1;

      wig_table_init(2 * 100, 9);
      wig_temp_init(2 * 100);

      double val3j = wig3jj(ll1, ll2, ll3, mm1, mm2, -mm3);

      wig_temp_free();
      wig_table_free();

      val3j *= phase * std::sqrt(ll3 + 1);

      // std::cout << ll1 << "," << ll2 << "," << ll3 << "," << mm1 << "," <<
      // mm2 << "," << mm3 << "->" << val3j << std::endl;

      return val3j;
    }
  }

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
  void HighestWeightReducedCoeffs(
      const std::vector<u3::gelfand>& gArr1,
      const std::vector<u3::gelfand>& gArr2,
      const std::vector<u3::gelfand>& gArr3,
      std::vector<std::tuple<int, int, int>>& tArr,
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& coeffs)
  {

    std::vector<std::pair<int, int>> indices =
        expand_to_uncoupled(gArr3[0], gArr1, gArr2);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> mat =
        E23E12(gArr1, gArr2, indices);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nullspace =
        null_space(mat);

    for (int i = 0; i < indices.size(); i++)
    {
      int i1 = indices[i].first - gArr1[indices[i].first].n12() +
               gArr1[indices[i].first].n11();
      int i2 = indices[i].second - gArr2[indices[i].second].n12() +
               gArr2[indices[i].second].n11();

      std::tuple t(i1, i2, 0);
      auto it = std::find(tArr.begin(), tArr.end(), t);

      if (it == tArr.end())
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
          double cgc = u3::CGC_SO3(jj1, jj2, jj3, mm1, mm2, mm3);

          coeffs.row((i1 * gArr2.size() + i2) * gArr3.size()) =
              nullspace.row(i) / cgc;

          tArr.push_back(t);
        }
      }
    }

    return;
  }

  // lower-weight reduced coefficients
  // / n13a n23a n33a ; n13b n23b n33b || n13c n23c n33c \
  // \   n12a  n22a       n12b  n22b   ||   n12c n22c    /
  void LowerWeightReducedCoeffs(
      const std::vector<u3::gelfand>& gArr1,
      const std::vector<u3::gelfand>& gArr2,
      const std::vector<u3::gelfand>& gArr3,
      std::vector<std::tuple<int, int, int>>& tArr,
      std::vector<int>& same_weight_gArr3,
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& coeffs)
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

      // find parent states and operators to get to lower weight states
      // (daughter states)
      std::vector<std::pair<int, int>> actOnParent;

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

      // matrix b is independent of outer multiplicity
      Eigen::MatrixXd b(actOnParent.size(), same_weight_gArr3.size());

      // the lowering matrix from the action of E32 and E31, which depends on
      // outer multiplicity
      int rmax = coeffs.cols();

      std::vector<Eigen::MatrixXd> E32E31(rmax);

      // set matrix elements to zero
      for (int r = 0; r < rmax; r++)
      {
        E32E31[r] = Eigen::MatrixXd(actOnParent.size(),
                                    uncoupled_same_weight_gArr3.size());
        E32E31[r].setZero();
      }

      // now build matrix b and E32E31
      for (int i = 0; i < actOnParent.size(); i++)
      {
        int opt = actOnParent[i].first, parent_g3 = actOnParent[i].second;

        int parent_g3_hw =
            parent_g3 - gArr3[parent_g3].n12() + gArr3[parent_g3].n11();

        std::vector<std::pair<int, int>> uncoupled_parent_g3 =
            expand_to_uncoupled(gArr3[parent_g3], gArr1, gArr2);

        // build with E32
        if (opt == 0)
        {
          // loop over daughter states
          for (int j = 0; j < same_weight_gArr3.size(); j++)
            if (u3::E23_12(gArr3[same_weight_gArr3[j]], gArr3[parent_g3]))
              b(i, j) = gArr3[parent_g3].e32_12();
            else if (u3::E23_22(gArr3[same_weight_gArr3[j]], gArr3[parent_g3]))
              b(i, j) = gArr3[parent_g3].e32_22();
            else
              b(i, j) = 0.0;

          for (int j = 0; j < uncoupled_same_weight_gArr3.size(); j++)
          {
            int g1 = uncoupled_same_weight_gArr3[j].first,
                g2 = uncoupled_same_weight_gArr3[j].second;

            for (int l = 0; l < uncoupled_parent_g3.size(); l++)
            {
              int parent_g1 = uncoupled_parent_g3[l].first,
                  parent_g2 = uncoupled_parent_g3[l].second;

              int parent_g1_hw = parent_g1 - gArr1[parent_g1].n12() +
                                 gArr1[parent_g1].n11(),
                  parent_g2_hw = parent_g2 - gArr2[parent_g2].n12() +
                                 gArr2[parent_g2].n11();

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
                int jj1 = gArr1[parent_g1].Lambdax2(),
                    jj2 = gArr2[parent_g2].Lambdax2(),
                    jj3 = gArr3[parent_g3].Lambdax2(),
                    mm1 = gArr1[parent_g1].MLambdax2(),
                    mm2 = gArr2[parent_g2].MLambdax2(),
                    mm3 = gArr3[parent_g3].MLambdax2();

                if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
                    (jj1 + jj2) >= jj3)
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

                  double cgc = u3::CGC_SO3(jj1, jj2, jj3, mm1, mm2, mm3);

                  for (int r = 0; r < rmax; r++)
                    E32E31[r](i, j) +=
                        coeffs((parent_g1_hw * gArr2.size() + parent_g2_hw) *
                                       gArr3.size() +
                                   parent_g3_hw,
                               r) *
                        cgc * factor;
                }
              }
            }
          }
        }

        // build with E31
        else if (opt == 1)
        {
          // loop over daughter states
          for (int j = 0; j < same_weight_gArr3.size(); j++)
            if (u3::E13_12(gArr3[same_weight_gArr3[j]], gArr3[parent_g3]))
              b(i, j) = gArr3[parent_g3].e31_12();
            else if (u3::E13_22(gArr3[same_weight_gArr3[j]], gArr3[parent_g3]))
              b(i, j) = gArr3[parent_g3].e31_22();
            else
              b(i, j) = 0.0;

          for (int j = 0; j < uncoupled_same_weight_gArr3.size(); j++)
          {
            int g1 = uncoupled_same_weight_gArr3[j].first,
                g2 = uncoupled_same_weight_gArr3[j].second;

            for (int l = 0; l < uncoupled_parent_g3.size(); l++)
            {
              int parent_g1 = uncoupled_parent_g3[l].first,
                  parent_g2 = uncoupled_parent_g3[l].second;

              int parent_g1_hw = parent_g1 - gArr1[parent_g1].n12() +
                                 gArr1[parent_g1].n11(),
                  parent_g2_hw = parent_g2 - gArr2[parent_g2].n12() +
                                 gArr2[parent_g2].n11();

              bool cond1 =
                  (u3::E13_12(gArr1[g1], gArr1[parent_g1]) && g2 == parent_g2);
              bool cond2 =
                  (u3::E13_22(gArr1[g1], gArr1[parent_g1]) && g2 == parent_g2);
              bool cond3 =
                  (u3::E13_12(gArr2[g2], gArr2[parent_g2]) && g1 == parent_g1);
              bool cond4 =
                  (u3::E13_22(gArr2[g2], gArr2[parent_g2]) && g1 == parent_g1);

              if (cond1 || cond2 || cond3 || cond4)
              {
                int jj1 = gArr1[parent_g1].Lambdax2(),
                    jj2 = gArr2[parent_g2].Lambdax2(),
                    jj3 = gArr3[parent_g3].Lambdax2(),
                    mm1 = gArr1[parent_g1].MLambdax2(),
                    mm2 = gArr2[parent_g2].MLambdax2(),
                    mm3 = gArr3[parent_g3].MLambdax2();

                if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
                    (jj1 + jj2) >= jj3)
                {
                  double factor;
                  if (cond1)
                    factor = gArr1[parent_g1].e31_12();
                  else if (cond2)
                    factor = gArr1[parent_g1].e31_22();
                  else if (cond3)
                    factor = gArr2[parent_g2].e31_12();
                  else
                    factor = gArr2[parent_g2].e31_22();

                  double cgc = u3::CGC_SO3(jj1, jj2, jj3, mm1, mm2, mm3);

                  for (int r = 0; r < rmax; r++)
                    E32E31[r](i, j) +=
                        coeffs((parent_g1_hw * gArr2.size() + parent_g2_hw) *
                                       gArr3.size() +
                                   parent_g3_hw,
                               r) *
                        cgc * factor;
                }
              }
            }
          }
        }
      }

      Eigen::ColPivHouseholderQR<Eigen::MatrixXd> factorize;
      factorize.compute(b);

      std::vector<Eigen::MatrixXd> coeffs_temp(rmax);

      for (int r = 0; r < rmax; r++)
        coeffs_temp[r] = factorize.solve(E32E31[r]);

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

          std::tuple t(i1_hw, i2_hw, i3_hw);
          auto it = std::find(tArr.begin(), tArr.end(), t);

          if (it == tArr.end())
          {
            int jj1 = gArr1[i1].Lambdax2(), jj2 = gArr2[i2].Lambdax2(),
                jj3 = gArr3[i3].Lambdax2(), mm1 = gArr1[i1].MLambdax2(),
                mm2 = gArr2[i2].MLambdax2(), mm3 = gArr3[i3].MLambdax2();

            if ((mm1 + mm2) == mm3 && std::fabs(jj1 - jj2) <= jj3 &&
                (jj1 + jj2) >= jj3)
            {
              double cgc = u3::CGC_SO3(jj1, jj2, jj3, mm1, mm2, mm3);

              tArr.push_back(t);
              for (int r = 0; r < rmax; r++)
              {
                if (std::fabs(coeffs_temp[r](i, j)) < 1e-10)
                  coeffs_temp[r](i, j) = 0.0;

                coeffs((i1_hw * gArr2.size() + i2_hw) * gArr3.size() + i3_hw,
                       r) = coeffs_temp[r](i, j) / cgc;
              }
            }
          }
        }
      }
    }

    return;
  }

  // all coupling coefficients
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
  WigU3U2(const u3::irrep& u3A, const u3::irrep& u3B, u3::irrep& u3C,
          const int& rhomax, std::vector<u3::gelfand>& gArrA,
          std::vector<u3::gelfand>& gArrB, std::vector<u3::gelfand>& gArrC,
          std::vector<std::tuple<int, int, int>>& tArr)
  {
    gArrA = u3::gelfandBasis(u3A);
    gArrB = u3::gelfandBasis(u3B);
    gArrC = u3::gelfandBasis(u3C);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffs(
        gArrA.size() * gArrB.size() * gArrC.size(), rhomax);
    coeffs.setZero();

    int w2max = u3C.n13() + u3C.n23();
    int w2min = u3C.n23() + u3C.n33();

    std::vector<std::vector<int>> all_same_weight_gArrC(w2max - w2min + 1);
    for (int i = 0; i < gArrC.size(); i++)
      if (gArrC[i].n12() == gArrC[i].n11())
        all_same_weight_gArrC[w2max - gArrC[i].n12() - gArrC[i].n22()]
            .push_back(i);

    u3::HighestWeightReducedCoeffs(gArrA, gArrB, gArrC, tArr, coeffs);

    for (int w2 = w2max - 1; w2 >= w2min; w2--)
    {
      u3::LowerWeightReducedCoeffs(gArrA, gArrB, gArrC, tArr,
                                   all_same_weight_gArrC[w2max - w2], coeffs);
    }

    u3C.normalize();
    for (int i = 0; i < gArrC.size(); i++)
      gArrC[i].normalize();

    return coeffs;
  }

} // namespace u3

#endif
