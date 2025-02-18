#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <decomposition.h>
#include <U3basis.h>
#include <U3U2xU3U2.h>

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

  auto [u3cs, rhomaxs] = u3::decomposition(u3a, u3b);

  std::cout << "The normalized coupled irreps are: " << std::endl;
  for (int i = 0; i < u3cs.size(); i++)
  {
    std::cout << "(" << i << "): [" << u3cs[i].n13() - u3cs[i].n33() << ","
              << u3cs[i].n23() - u3cs[i].n33() << ",0] \t rho = 1 ... "
              << rhomaxs[i] << std::endl;
  }

  int index;
  std::cout
      << "Please enter the index of the coupled irrep you are interested: ";
  std::cin >> index;

  int n33c = u3cs[index].n33();

  std::vector<u3::gelfand> gArrA, gArrB, gArrC;

  std::vector<std::tuple<int, int, int>> tArr;

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffs = u3::WigU3U2(
      u3a, u3b, u3cs[index], rhomaxs[index], gArrA, gArrB, gArrC, tArr);
  // after this, u3cs[index] is normalized, all gelfand states in gArrC are also
  // normalized

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < tArr.size(); i++)
  {
    int i1 = std::get<0>(tArr[i]), i2 = std::get<1>(tArr[i]),
        i3 = std::get<2>(tArr[i]);

    std::cout << "< " << gArrA[i1].n13() << "," << gArrA[i1].n23() << ",0;"
              << gArrA[i1].n12() << "," << gArrA[i1].n22() << " ; "
              << gArrB[i2].n13() << "," << gArrB[i2].n23() << ",0;"
              << gArrB[i2].n12() << "," << gArrB[i2].n22() << " || "
              << gArrC[i3].n13() << "," << gArrC[i3].n23() << ",0;"
              << gArrC[i3].n12() << "," << gArrC[i3].n22() << " >" << "\t"
              << coeffs.row((i1 * gArrB.size() + i2) * gArrC.size() + i3)
              << std::endl;
  }

  // Calculate the duration
  std::chrono::duration<double> duration = t1 - t0;

  // Output the time in seconds
  std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;

  return 0;
}
