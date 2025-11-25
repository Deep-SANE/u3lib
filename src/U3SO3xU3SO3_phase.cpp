#include <U3decomposition.h>
#include <U3basis.h>
#include <U3SO3xU3SO3.h> 

#ifdef USE_SLEPC
#include <slepcsys.h>
#endif 

auto diff(int k1, int k2, int k3,  
  const Eigen::MatrixXd& A, const Eigen::MatrixXd& B)
{
  Eigen::MatrixXd d(A.rows(), A.cols());

  for (size_t i = 0; i < k1; i++)
  {
  	for (size_t j = 0; j < k2; j++)
  	{
      for (size_t k = 0; k < k3; k++)
      {
        auto A_row = (i * k2 + j) * k3 + k;
        auto B_row = (k * k2 + j) * k1 + i;

        d.row(A_row) = A.row(A_row) - B.row(B_row);
      }
  	}
  }

  return d;
}

int main(int argc, char* argv[])
{
#ifdef USE_SLEPC
  SlepcInitialize(&argc, &argv, NULL, NULL);
#endif

#ifdef ADD_SU3LIB_PHASE_U3SO3

  qcx_utils::SO3_CGC cgcs;
  u3::WigU3U2_db wigU3U2db(cgcs);
  u3::WigU3SO3_db wigU3SO3db(cgcs, wigU3U2db);

  int lam1, mu1, lam2, mu2;
  std::cout << "Please enter lambda1, mu1: ";
  std::cin >> lam1 >> mu1;
  
  std::cout << "Please enter lambda2 != mu2: ";
  std::cin >> lam2 >> mu2;

  u3::irrep u3A(lam1, mu1), u3B(lam2, mu2);

  auto [u3Cs, rhomaxs] = u3::decomposition(u3A, u3B);

  std::cout << "The coupled irreps (lambda3,mu3) are" << std::endl;
  for (int i = 0; i < u3Cs.size(); i++)
  {
    std::cout << "(" << i << "): (" << u3Cs[i].lam() << "," << u3Cs[i].mu()
              << ") \t rho = 1 ... " << rhomaxs[i] << std::endl;
  }

  int index;
  std::cout << "Please enter the index of the irrep you are interested: ";
  std::cin >> index;
  if (rhomaxs[index] > 1)
  {
    std::cout << "Only rhomax = 1 can be tested";
    return 1;
  }

  u3::irrep u3C = u3Cs[index];
  int lam3 = u3C.lam(), mu3 = u3C.mu(), rhomax = rhomaxs[index];

  auto pArrA = u3::physicalBasis(u3A), 
  	   pArrB = u3::physicalBasis(u3B), 
       pArrC = u3::physicalBasis(u3C);

  double dimA = u3A.dim(), dimC = u3C.dim();

  for (const auto& pA : pArrA)
  {
  	for (const auto& pB : pArrB)
  	{
  	  for (const auto& pC : pArrC)
  	  {
  	  	auto [flag, cgcU3] = wigU3SO3db.get(lam1, mu1, pA.L,
  	  										   lam2, mu2, pB.L, lam3, mu3, pC.L);

  	  	if (flag)
  	  	{
  	  	  auto [flag_conj, cgcU3_conj] = wigU3SO3db.get(lam3, mu3, pC.L,
  	  										               mu2, lam2, pB.L, lam1, mu1, pA.L);
          assert(flag_conj);

          int phi = lam1 + mu1 - pA.L - (lam3 + mu3 - pC.L) - pB.L;

  	  	  double phase = (phi % 2 == 0) ? 1.0 : -1.0;
  	  	  double factor = 
  	  	  	std::sqrt((2.0 * pA.L + 1.0) * dimC / (2.0 * pC.L + 1.0) / dimA);
          
          auto d = diff(pA.kappa, pB.kappa, pC.kappa, cgcU3, phase * factor * cgcU3_conj);

  	  	  std::cout << "\nl1 = " << pA.L << "; l2 = " << pB.L << "; l3 = " << pC.L 
  	  	  			<< "\nCoupling:\n" << cgcU3 
  	  	  			<< "\nConjugate coupling:\n" << cgcU3_conj 
  	  	  			<< "\nDifference with factor: " << ((d.norm() > 1e-10) ? "Violated" : "")
                //<< "\nDifference with factor: " << ((d.norm() > 1e-6) ? d : Eigen::MatrixXd(0,0))
  	  	  			<< std::endl;
  	  	}
  	  }
  	}
  }

#else
  std::cout << "No phase enforced" << std::endl;
#endif

#ifdef USE_SLEPC 
  SlepcFinalize();
#endif

  return 0;
}
