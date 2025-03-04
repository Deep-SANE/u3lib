#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <decomposition.h>
#include <U3basis.h>
#include <U3SO3xU3SO3_specific.h>

int main()
{

  int lam1, mu1, lam2, mu2;
  std::cout << "Please enter lambda1, mu1: ";
  std::cin >> lam1 >> mu1;

  std::cout << "Please enter lambda2, mu2: "; 
  std::cin >> lam2 >> mu2;

  u3::irrep u3a(lam1, mu1), u3b(lam2, mu2);

  auto [u3cs, rhomaxs] = u3::decomposition(u3a, u3b);

  std::cout << "The coupled irreps (lambda,mu) are" << std::endl;
  for (int i = 0; i < u3cs.size(); i++)
  {
    std::cout << "(" << i << "): (" << u3cs[i].lam() << "," << u3cs[i].mu()
              << ") \t rho = 1 ... " << rhomaxs[i] << std::endl;
  }

  int index;
  std::cout << "Please enter the index of the irrep you are interested: ";
  std::cin >> index;

  int l1, l2, l3;
  std::cout << "Please enter l1 <= (lambda1 + mu1): ";
  std::cin >> l1;
  std::cout << "Please enter l2 <= (lambda2 + mu2): ";
  std::cin >> l2;
  std::cout << "Please enter l3 <= (lambda3 + mu3): ";
  std::cin >> l3;

  if (std::fabs(l1 - l2) <= l3 && (l1 + l2) >= l3)
  {
    u3::physicalvector pA(u3a,l1), pB(u3b,l2), pC(u3cs[index],l3);

    // Start time measurement
    auto t0 = std::chrono::high_resolution_clock::now();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffs =
                  u3::WigU3SO3_specific(rhomaxs[index], pA, pB, pC);

    // End time measurement
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Reduced coupling coefficients: < (" << lam1 << "," << mu1 << ") k1 " << l1 
              << " ; (" << lam2 << "," << mu2 << ") k2 " << l2 << " || "
              << "(" << pC.toprow.lam() << "," << pC.toprow.mu() << ") k3 " << l3 << " > "
              << "for rho = 1 ... " << rhomaxs[index] << std::endl;
    std::cout << "k1\tk2\tk3" << std::endl;

    for (int k1 = 0; k1 < pA.kappa; k1++)
      for (int k2 = 0; k2 < pB.kappa; k2++)
        for (int k3 = 0; k3 < pC.kappa; k3++)
          std::cout << "" << k1+1 << "\t" << k2+1 << "\t" << k3+1 << "\t" 
              << coeffs.row((k1 * pB.kappa + k2) * pC.kappa + k3) << std::endl;  

    // Calculate the duration
    std::chrono::duration<double> duration = t1 - t0;

    // Output the time in seconds
    std::cout << "Elapsed time: " << duration.count() << " seconds or "
              << duration.count() / 60.0 << " minutes" << std::endl;
  }
  else
    std::cout << "Coupling does not exist" << std::endl;
  
  return 0;
}
 