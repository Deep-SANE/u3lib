#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include <U3basis.h>
#include <U3SO3_to_U3U2U1.h> 

#ifdef USE_SLEPC
#include <slepcsys.h>
#endif
  
int main(int argc, char* argv[])
{
#ifdef USE_SLEPC
  SlepcInitialize(&argc, &argv, NULL, NULL);
#endif

  int lambda, mu;
  std::cout << "Enter SU(3) quantum numbers (lambda,mu): ";
  std::cin >> lambda >> mu;

  u3::irrep x(lambda, mu);

  // building the gelfand basis
  std::vector<u3::gelfand> gArr = u3::gelfandBasis(x);

  // building the physical basis excluding those which have zero multiplicity
  std::vector<u3::physicalvector> pArr = u3::physicalBasis(x);

  std::cout << "The physical basis states are: " << std::endl;
  for (int i = 0; i < pArr.size(); i++)
  {
    std::cout << "(" << i << "): |" << "(" << pArr[i].lam() << ","
              << pArr[i].mu() << "), kappa_max=" << pArr[i].kappa
              << ", L=" << pArr[i].L << ">" << std::endl;
  } 

  int index;
  std::cout << "Enter the index of the state that you are interested: ";
  std::cin >> index;
  const auto& p = pArr[index];

  int mode;
  std::cout << "Enter 0 for specific projection M, 1 for all M's: ";
  std::cin >> mode;

  if (mode == 0)
  {
    int M;
    std::cout << "Enter projection M <= L: ";
    std::cin >> M;

    // Start time measurement
    auto t0 = std::chrono::high_resolution_clock::now();

    // get the transformation
    auto [canoBasis, transformation] = u3::phys_to_cano(p, M, gArr);

    // End time measurement
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout
      << "The Gelfand states and the expansion coefficients for kappa = 1,...,"
      << p.kappa << " are:" << std::endl;

    for (int i = 0; i < canoBasis.size(); i++)
      std::cout << "|" << gArr[canoBasis[i]] << "> \t" << transformation.row(i)
              << std::endl;

    // Calculate the duration
    std::chrono::duration<double> duration = t1 - t0;

    // Output the time in seconds
    std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;
  }
  else if (mode == 1)
  {
    // Start time measurement
    auto t0 = std::chrono::high_resolution_clock::now();

    // get the transformation
    auto [canoBasis, transformation] = u3::phys_to_cano_allM(p, gArr);

    // End time measurement
    auto t1 = std::chrono::high_resolution_clock::now();

    for (int m = p.L; m >= -p.L; m--)
    {
      std::cout
        << "Gelfand states and expansion coefficients for M = " << m 
        << " for kappa = 1,...," << p.kappa << " are:" << std::endl;

      for (int i = 0; i < canoBasis[p.L - m].size(); i++)
        std::cout << "|" << gArr[canoBasis[p.L - m][i]] 
                  << "> \t" << transformation[p.L - m].row(i)
                  << std::endl;
    }

    // Calculate the duration
    std::chrono::duration<double> duration = t1 - t0;

    // Output the time in seconds
    std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;
  }
  else
  {
    std::cout << "Wrong mode" << std::endl;
  }

#ifdef USE_SLEPC
  SlepcFinalize();
#endif

  return 0;
}
