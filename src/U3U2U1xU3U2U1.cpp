#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <decomposition.h>
#include <U3basis.h>
#include <U3U2U1xU3U2U1.h>

int main()
{

  int n13a, n23a, n33a;
  std::cout << "Please enter the first irrep [n13a,n23a,n33a]: ";
  std::cin >> n13a >> n23a >> n33a;
  u3::irrep u3a(n13a, n23a, n33a);
  if (n33a != 0)
  {
    u3a.normalize();
    std::cout << "Your irrep is normalized to " << u3a << std::endl;
  }

  int n13b, n23b, n33b;
  std::cout << "Please enter the second irrep [n13b,n23b,n33b]: ";
  std::cin >> n13b >> n23b >> n33b;
  u3::irrep u3b(n13b, n23b, n33b);
  if (n33b != 0)
  {
    u3b.normalize();
    std::cout << "Your irrep is normalized to " << u3b << std::endl;
  }

  std::vector<int> rhomaxs;
  std::vector<u3::irrep> u3cs = u3::decomposition(u3a, u3b, rhomaxs);

  std::cout << "The normalized coupled irreps are: " << std::endl;
  for (int i = 0; i < u3cs.size(); i++)
  {
    std::cout << "(" << i << "): [" << u3cs[i].n13 - u3cs[i].n33 << ","
              << u3cs[i].n23 - u3cs[i].n33 << ",0] \t rho = 1 ... "
              << rhomaxs[i] << std::endl;
  }

  int index;
  std::cout
      << "Please enter the index of the coupled irrep you are interested: ";
  std::cin >> index;

  int n33c = u3cs[index].n33;

  std::vector<u3::gelfand> gArrA, gArrB, gArrC;

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffs =
      u3::WigU3U2U1(u3a, u3b, u3cs[index], rhomaxs[index], gArrA, gArrB, gArrC);
  // after this, u3cs[index] is normalized, all gelfand states in gArrC are also
  // normalized

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < gArrA.size(); i++)
  {
    for (int j = 0; j < gArrB.size(); j++)
    {
      for (int k = 0; k < gArrC.size(); k++)
      {
        // print out only coupling coefficients allowed by selection rule
        if ((gArrA[i].e11() + gArrB[j].e11()) == (gArrC[k].e11() + n33c) &&
            (gArrA[i].e22() + gArrB[j].e22()) == (gArrC[k].e22() + n33c))
          std::cout << "< " << gArrA[i] << " ; " << gArrB[j] << " | "
                    << gArrC[k] << " >" << "\t"
                    << coeffs.row((i * gArrB.size() + j) * gArrC.size() + k)
                    << std::endl;
      }
    }
  }

  // Calculate the duration
  std::chrono::duration<double> duration = t1 - t0;

  // Output the time in seconds
  std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;

  return 0;
}
