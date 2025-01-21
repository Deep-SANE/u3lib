#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include <U3basis.h>
#include <U3SO3_to_U3U2U1.h>
#include <utils/nullspaceSolver.h>

int main()
{

  int lambda, mu;
  std::cout << "Please enter SU(3) quantum numbers (lambda,mu): ";
  std::cin >> lambda >> mu;

  u3::irrep u3top(lambda, mu);

  // building the gelfand basis
  std::vector<u3::gelfand> gArr = u3::gelfandBasis(u3top);

  // building the physical basis excluding those which have zero multiplicity
  std::vector<u3::physicalvector> pArr = u3::physicalBasis(u3top);

  std::cout << "The physical basis states are: " << std::endl;
  for (int i = 0; i < pArr.size(); i++)
  {
    std::cout << "(" << i << "): |" << "(" << pArr[i].toprow.lam() << ","
              << pArr[i].toprow.mu()
              << "), kappa_max=" << pArr[i].kappa
              << ", L=" << pArr[i].L << ">" << std::endl;
  }

  int index;
  std::cout << "Please enter the index of the state that you are interested: ";
  std::cin >> index;

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  // store Gelfand states associated with the interested physical state
  std::vector<int> indices;

  // get the transformation
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nullspace = u3::phys_to_cano(pArr[index], indices, gArr);

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout
      << "The Gelfand states and the expansion coefficients for kappa = 1,...,"
      << pArr[index].kappa << " are:" << std::endl;

  for (int i = 0; i < indices.size(); i++)
    std::cout << "|" << gArr[indices[i]] << "> \t" << nullspace.row(i)
              << std::endl;

  // Calculate the duration
  std::chrono::duration<double> duration = t1 - t0;

  // Output the time in seconds
  std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;

  // Print error
  // std::cout << "Error: " << error(Lproj,nullspace) << std::endl; 

  // now going lower from the HW state, remember that gamma = L - M = 0,1,...,2L
  // store transformation, each column corresponds to one kappa
      std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
          allCanoCoeffs(2*pArr[index].L+1);
      // store expansion
      std::vector<std::vector<int>> allCanoBasis(2*pArr[index].L+1); 

      allCanoCoeffs[0] = nullspace;
      allCanoBasis[0] = indices;

      for (int gamma = 1; gamma <= 2 * pArr[index].L; gamma++)
      {

        // store new expansion
        allCanoBasis[gamma] = u3::expand_to_gelfand(pArr[index].L - gamma, gArr);

        // normalize coefficient
        double n = std::sqrt(2.0 / (gamma * (2 * pArr[index].L - gamma + 1.0)));

        allCanoCoeffs[gamma] = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>(allCanoBasis[gamma].size(), pArr[index].kappa);
        allCanoCoeffs[gamma].setZero();

        // basis states are obtained by acting with two generators E31 and E23
        for (int m = 0; m < allCanoBasis[gamma].size(); m++)
        {
          for (int j=0; j<allCanoBasis[gamma-1].size(); j++)
          {
            if ( u3::E13_12(gArr[allCanoBasis[gamma][m]],gArr[allCanoBasis[gamma-1][j]]) ){
              for (int k = 0; k < pArr[index].kappa; k++)
                allCanoCoeffs[gamma](m,k) += n*allCanoCoeffs[gamma-1](j,k)*gArr[allCanoBasis[gamma-1][j]].e31_12();
            }
            else if ( u3::E13_22(gArr[allCanoBasis[gamma][m]],gArr[allCanoBasis[gamma-1][j]]) ){
              for (int k = 0; k < pArr[index].kappa; k++)
                allCanoCoeffs[gamma](m,k) += n*allCanoCoeffs[gamma-1](j,k)*gArr[allCanoBasis[gamma-1][j]].e31_22();
            }
            else if ( u3::E23_12(gArr[allCanoBasis[gamma-1][j]],gArr[allCanoBasis[gamma][m]]) ){
              for (int k = 0; k < pArr[index].kappa; k++)
                allCanoCoeffs[gamma](m,k) += n*allCanoCoeffs[gamma-1](j,k)*gArr[allCanoBasis[gamma-1][j]].e23_12();
            }
            else if ( u3::E23_22(gArr[allCanoBasis[gamma-1][j]],gArr[allCanoBasis[gamma][m]]) ){
              for (int k = 0; k < pArr[index].kappa; k++)
                allCanoCoeffs[gamma](m,k) += n*allCanoCoeffs[gamma-1](j,k)*gArr[allCanoBasis[gamma-1][j]].e23_22();
            }
          }
        }

        std::cout
      << "The Gelfand states and the expansion coefficients for M = " << pArr[index].L - gamma << ", kappa = 1,...,"
      << pArr[index].kappa << " are:" << std::endl;

      for (int i = 0; i < allCanoBasis[gamma].size(); i++)
          std::cout << "|" << gArr[allCanoBasis[gamma][i]] << "> \t" << allCanoCoeffs[gamma].row(i)
              << std::endl;
      }

  return 0;
}
