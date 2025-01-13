#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <decomposition.h>
#include <U3basis.h>
#include <U3SO3xU3SO3.h>

int main()
{

  int lam1, mu1, lam2, mu2;
  std::cout << "Please enter lambda1, mu1: ";
  std::cin >> lam1 >> mu1;

  std::cout << "Please enter lambda2, mu2: ";
  std::cin >> lam2 >> mu2;

  u3::irrep u3a(lam1, mu1), u3b(lam2, mu2);

  std::vector<int> rhomaxs;
  std::vector<u3::irrep> u3cs = u3::decomposition(u3a, u3b, rhomaxs);

  std::cout << "The coupled irreps (lambda,mu) are" << std::endl;
  for (int i = 0; i < u3cs.size(); i++)
  {
    std::cout << "(" << i << "): (" << u3cs[i].lam() << "," << u3cs[i].mu()
              << ") \t rho = 1 ... " << rhomaxs[i] << std::endl;
  }

  int index;
  std::cout << "Please enter the index of the irrep you are interested: ";
  std::cin >> index;

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  std::vector<u3::physicalvector> physvecA;
  std::vector<u3::physicalvector> physvecB;
  std::vector<u3::physicalvector> physvecC;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffs =
      u3::WigU3SO3(u3a, u3b, u3cs[index], physvecA, physvecB, physvecC);

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();

  for (int j = 0; j < physvecA.size(); j++)
  {
    std::cout << "< " << physvecA[j] << " ; " << physvecB[j] << " || "
              << physvecC[j] << " >" << "\t" << coeffs.row(j) << std::endl;
  }

  // Calculate the duration
  std::chrono::duration<double> duration = t1 - t0;

  // Output the time in seconds
  std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;

  return 0;
}
