#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple> 
#include <vector> 

#include <U3decomposition.h>
#include <U3basis.h>
#include <U3SO3xU3SO3.h> 
 
#ifdef USE_SLEPC
#include <slepcsys.h>
#endif 

int main(int argc, char* argv[])
{
#ifdef USE_SLEPC
  SlepcInitialize(&argc, &argv, NULL, NULL);
#endif

  int lam1, mu1, lam2, mu2;
  std::cout << "Please enter lambda1, mu1: ";
  std::cin >> lam1 >> mu1;

  std::cout << "Please enter lambda2, mu2: ";
  std::cin >> lam2 >> mu2;

  u3::irrep u3A(lam1, mu1), u3B(lam2, mu2);

  auto [u3Cs, rhomaxs] = u3::decomposition(u3A, u3B);

  std::cout << "The coupled irreps (lambda,mu) are" << std::endl;
  for (int i = 0; i < u3Cs.size(); i++)
  {
    std::cout << "(" << i << "): (" << u3Cs[i].lam() << "," << u3Cs[i].mu()
              << ") \t rho = 1 ... " << rhomaxs[i] << std::endl;
  }

  int index;
  std::cout << "Please enter the index of the irrep you are interested: ";
  std::cin >> index;
  u3::irrep u3C = u3Cs[index];
  int rhomax = rhomaxs[index];

  // Build the physical bases
  std::vector<u3::physicalvector> pArrA, pArrB, pArrC;

  // Input angular momentum
  int l1, l2, l3;
  std::cout << "Please enter angular momentum values 0 <= l <= (lambda + mu), "
               "(enter -1 to stop, -2 to include all)"
            << std::endl;

  std::cout << "l1: ";
  while (true)
  {
    std::cin >> l1;

    // Stop condition
    if (l1 == -1) break;

    // Include all
    if (l1 == -2)
    {
      pArrA = u3::physicalBasis(u3A);
      break;
    }

    // Normal case: Add valid values to vector
    if (l1 >= 0)
    {
      u3::physicalvector p(u3A, l1);
      if (p.kappa() > 0) pArrA.push_back(p);
    }
  }

  std::cout << "l2: ";
  while (true)
  {
    std::cin >> l2;

    // Stop condition
    if (l2 == -1) break;

    // Include all
    if (l2 == -2)
    {
      pArrB = u3::physicalBasis(u3B);
      break;
    }

    // Normal case: Add valid values to vector
    if (l2 >= 0)
    {
      u3::physicalvector p(u3B, l2);
      if (p.kappa() > 0) pArrB.push_back(p);
    }
  }

  std::cout << "l3: ";
  while (true)
  {
    std::cin >> l3;

    // Stop condition
    if (l3 == -1) break;

    // Include all
    if (l3 == -2)
    {
      pArrC = u3::physicalBasis(u3C);
      break;
    }

    // Normal case: Add valid values to vector
    if (l3 >= 0)
    {
      u3::physicalvector p(u3C, l3);
      if (p.kappa() > 0) pArrC.push_back(p);
    }
  }
 
  if (pArrA.empty() || pArrB.empty() || pArrC.empty())
  {
    std::cout << "Your angular momentum values are invalid" << std::endl;
    return 1;
  }

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  qcx_utils::SO3_CGC cgcs;
  u3::WigU3U2_db wigU3U2db(cgcs);

  auto [tArr, coeffs] = 
    u3::WigU3SO3(rhomax, pArrA, pArrB, pArrC, cgcs, wigU3U2db);

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();

  // Print out
  for (int t = 0; t < tArr.size(); t++)
  {
    int i = std::get<0>(tArr[t]), j = std::get<1>(tArr[t]),
        k = std::get<2>(tArr[t]);

    std::cout << "Reduced coupling coefficients: < (" << lam1 << "," << mu1
              << ") k1 " << pArrA[i].L() << " ; (" << lam2 << "," << mu2
              << ") k2 " << pArrB[j].L() << " || "
              << "(" << u3C.lam() << "," << u3C.mu() << ") k3 " << pArrC[k].L()
              << " > "
              << "for rho = 1 ... " << rhomax << std::endl;
    std::cout << "k1\tk2\tk3" << std::endl;

    Eigen::MatrixXd coeffs_temp = coeffs[t];

    for (int k1 = 0; k1 < pArrA[i].kappa(); k1++)
      for (int k2 = 0; k2 < pArrB[j].kappa(); k2++)
        for (int k3 = 0; k3 < pArrC[k].kappa(); k3++)
          std::cout << "" << k1 + 1 << "\t" << k2 + 1 << "\t" << k3 + 1 << "\t"
                    << coeffs_temp.row(
                           (k1 * pArrB[j].kappa() + k2) * pArrC[k].kappa() + k3)
                    << std::endl;
  }

  // Calculate the duration
  std::chrono::duration<double> duration = t1 - t0;

  // Output the time in seconds
  std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;

#ifdef USE_SLEPC 
  SlepcFinalize();
#endif

  return 0;
}
