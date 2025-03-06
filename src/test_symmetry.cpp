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
  std::vector<u3::physicalvector> pArr = u3::physicalBasis(u3top);

/*
  std::cout << "The physical basis states are: " << std::endl;
  for (int i = 0; i < pArr.size(); i++)
  {
    std::cout << "(" << i << "): |" << "(" << pArr[i].toprow.lam() << ","
              << pArr[i].toprow.mu() << "), kappa_max=" << pArr[i].kappa
              << ", L=" << pArr[i].L << ">" << std::endl;
  }
*/

  for (int index = 0; index < pArr.size(); index++)
  {
    std::vector<std::vector<int>> canoBasis(2 * pArr[index].L + 1);
    std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> transformation(2 * pArr[index].L + 1);
    transformation[0] = u3::phys_to_cano(pArr[index], pArr[index].L, canoBasis[0], gArr);

    // ensure the nullspace is solved correctly
    assert(pArr[index].kappa == transformation[0].cols());

    for (int m = pArr[index].L - 1; m >= -pArr[index].L; m--)
    {
      canoBasis[pArr[index].L - m] = u3::expand_to_gelfand(m, gArr);
      transformation[pArr[index].L - m] = u3::lowerM(pArr[index].kappa, pArr[index].L, m, gArr, 
        transformation[pArr[index].L - m - 1], canoBasis[pArr[index].L - m -1], canoBasis[pArr[index].L - m]);
    }

    double err = 0.0; 
 
    for (int i = 0; i < canoBasis.size(); i++)
    { 
      //std::cout << "For M = " << pArr[index].L - i << " for kappa = 1,...," 
      //                          << pArr[index].kappa << " are:" << std::endl;

      if (i > pArr[index].L)
      {
        //for (int j = 0; j < canoBasis[i].size(); j++) 
        // std::cout << "|" << gArr[canoBasis[i][j]] << "> \t" << transformation[i].row(j) << std::endl;
      //}
      //else
      //{
        for (int j = 0; j < canoBasis[i].size(); j++) 
        {
          int phase = ((lambda + pArr[index].L + gArr[canoBasis[i][j]].n22())%2 == 0) ? 1 : -1;

          Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> diff = 
            transformation[i].row(j) -  phase*transformation[2 * pArr[index].L - i].row(j);

          for (int k = 0; k < diff.cols(); k++)
            err += std::fabs(diff(0,k));

          //std::cout << "|" << gArr[canoBasis[i][j]] << "> \t" << diff << std::endl;
        }
      }
    }
  
    // Print error
    std::cout << "Error: " << err << std::endl;
  }  
   
  return 0;
}