#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include <U3basis.h>
#include <U3SO3_to_U3U2U1.h> 

#ifdef U3LIB_USE_SLEPC
#include <slepcsys.h>
#endif
  
int main(int argc, char* argv[])
{
#ifdef U3LIB_USE_SLEPC
  SlepcInitialize(&argc, &argv, NULL, NULL);
#endif

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
              << pArr[i].toprow.mu() << "), kappa_max=" << pArr[i].kappa
              << ", L=" << pArr[i].L << ">" << std::endl;
  } 

  int index;
  std::cout << "Please enter the index of the state that you are interested: ";
  std::cin >> index;

  int M;
  std::cout << "Please enter projection M <= L: ";
  std::cin >> M;

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  // store Gelfand states associated with the interested physical state
  std::vector<int> canoBasis;

  // get the transformation
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> transformation =
      u3::phys_to_cano(pArr[index], M, canoBasis, gArr);

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout
      << "The Gelfand states and the expansion coefficients for kappa = 1,...,"
      << pArr[index].kappa << " are:" << std::endl;

  for (int i = 0; i < canoBasis.size(); i++)
    std::cout << "|" << gArr[canoBasis[i]] << "> \t" << transformation.row(i)
              << std::endl;

  // Calculate the duration
  std::chrono::duration<double> duration = t1 - t0;

  // Output the time in seconds
  std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;

  // Print error
  // std::cout << "Error: " << error(Lproj,nullspace) << std::endl;

#ifdef U3LIB_USE_SLEPC
  SlepcFinalize();
#endif

  return 0;
}
