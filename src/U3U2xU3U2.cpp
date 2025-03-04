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

  std::vector<u3::gelfand> gArrA, gArrB, gArrC;

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  auto coeffs =
      u3::WigU3U2(u3a, u3b, u3cs[index], rhomaxs[index], gArrA, gArrB, gArrC);
  // after this, u3cs[index] is normalized, all gelfand states in gArrC are also
  // normalized

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();

  const auto& u3c = u3cs[index];

  for (int n12a = u3a.n13(); n12a >= u3a.n23(); n12a--)
    for (int n22a = u3a.n23(); n22a >= u3a.n33(); n22a--)
      for (int n12b = u3b.n13(); n12b >= u3b.n23(); n12b--)
        for (int n22b = u3b.n23(); n22b >= u3b.n33(); n22b--)
          for (int n12c = u3c.n13(); n12c >= u3c.n23(); n12c--)
            for (int n22c = u3c.n23(); n22c >= u3c.n33(); n22c--)
            {
              auto acc = coeffs.access(n12a, n22a, n12b, n22b, n12c, n22c);
              
              if (acc.allowed == false)
                continue;

              std::cout << "< " << u3a.n13() << "," << u3a.n23() << ",0;"
                        << n12a << "," << n22a << " ; " << u3b.n13() << ","
                        << u3b.n23() << ",0;" << n12b << "," << n22b << " || "
                        << u3c.n13() << "," << u3c.n23() << ",0;" << n12c << ","
                        << n22c << " >" << "\t" << coeffs.get_row(acc)
                        << std::endl;
            } 
 
  std::cout << "Reduced coefficients count: " << coeffs.size() << std::endl;

  // Calculate the duration
  std::chrono::duration<double> duration = t1 - t0;

  // Output the time in seconds
  std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;

  return 0;
}
