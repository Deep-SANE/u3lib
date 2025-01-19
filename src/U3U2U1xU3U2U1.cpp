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
  std::vector<u3::gelfand> gArr1 = gelfandBasis(u3a);

  int n13b, n23b, n33b;
  std::cout << "Please enter the second irrep [n13b,n23b,n33b]: ";
  std::cin >> n13b >> n23b >> n33b;
  u3::irrep u3b(n13b, n23b, n33b);
  if (n33b != 0)
  {
    u3b.normalize();
    std::cout << "Your irrep is normalized to " << u3b << std::endl;
  }
  std::vector<u3::gelfand> gArr2 = gelfandBasis(u3b);

  std::vector<int> rhomaxs;
  std::vector<u3::irrep> u3cs = decomposition(u3a, u3b, rhomaxs);

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

  std::vector<u3::gelfand> gArr3 = gelfandBasis(u3cs[index]);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> coeffs(
      u3a.dim() * u3b.dim() * u3cs[index].dim(), rhomaxs[index]);
  for (int i = 0; i < coeffs.rows(); i++)
    for (int j = 0; j < coeffs.cols(); j++)
      coeffs(i, j) = 0.0;

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  for (int w3 = u3cs[index].n33; w3 <= u3cs[index].n13; w3++)
  {
    for (int w1 = u3cs[index].n13; w1 >= u3cs[index].n33; w1--)
    {
      if (w3 == u3cs[index].n33 && w1 == u3cs[index].n13)
      {
      u3:
        HighestWeightCoeffs(gArr1, gArr2, gArr3, coeffs);
      }
      else
      {
        u3::LowerWeightCoeffs(w3, w1, gArr1, gArr2, gArr3, coeffs);
      }
    }
  }

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();
  /*
    for (int i=0; i<gArr1.size(); i++)
    {
      for (int j=0; j<gArr2.size(); j++)
      {
        for (int k=0; k<gArr3.size(); k++)
        {
          if (i==0 && j==0) { gArr3[k].normalize(); }
          std::cout << "< " << gArr1[i] << " ; " << gArr2[j] << " | " <<
    gArr3[k]
                << " >" << "\t" << coeffs.row( (i*gArr2.size()+j)*gArr3.size()+k
    ) << std::endl;
        }
      }
    }
  */
  // Calculate the duration
  std::chrono::duration<double> duration = t1 - t0;

  // Output the time in seconds
  std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;

  return 0;
}
