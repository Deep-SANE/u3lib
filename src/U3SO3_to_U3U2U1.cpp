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
  std::vector<u3::physicalvector> physbasis = u3::physicalBasis(u3top);

  std::cout << "The physical basis states are: " << std::endl;
  for (int i = 0; i < physbasis.size(); i++)
  {
    std::cout << "(" << i << "): |" << "(" << physbasis[i].toprow.lam() << ","
              << physbasis[i].toprow.mu()
              << "), kappa_max=" << physbasis[i].kappa
              << ", L=" << physbasis[i].L << ">" << std::endl;
  }

  int index;
  std::cout << "Please enter the index of the state that you are interested: ";
  std::cin >> index;

  // find Gelfand states associated with the interested physical state
  std::vector<int> indices = u3::expand_to_gelfand(physbasis[index].L, gArr);

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  // Building the projection matrix
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Lproj =
      Lprojection(physbasis[index], gArr, indices);
  // and solving its null space
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> nullspace =
      phys_to_cano(physbasis[index], Lproj);
  // after solving the nullspace, the dimension of the nullspace is the
  // multiplicity of L, and the multiplicity within the physical basis state
  // labeling is modified

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout
      << "The Gelfand states and the expansion coefficients for kappa = 1,...,"
      << physbasis[index].kappa << " are:" << std::endl;

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

  return 0;
}
