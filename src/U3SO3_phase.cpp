#include <U3decomposition.h>
#include <U3basis.h>
#include <U3SO3xU3SO3.h> 

#ifdef USE_SLEPC
#include <slepcsys.h>
#endif 


/* This program tests the following conjugation relation
   
   < (lam1,mu1) k1 l1; (n,0) l || (lam3,mu3) k3 l3 >_rhomax=1
   = (-)^phi * sqrt(dim(lam3,mu3) * dim(l1) / dim(lam1,mu1) / dim(l3))
   * < (lam3,mu3) k3 l3; (0,n) l || (lam1,mu1) k1 l1 >_rhomax=1

   where
   phi = lam1 + mu1 - l1 - (lam3 + mu3 - l3) + (n + l)/2
*/

Eigen::MatrixXd diff(int k1, const Eigen::MatrixXd& m1, 
  					 int k3, const Eigen::MatrixXd& m3)
{
  Eigen::MatrixXd d(m1.rows(), m1.cols());

  for (size_t i = 0; i < k1; i++)
  {
  	for (size_t j = 0; j < k3; j++)
  	{
  	  auto m1_row = i * k3 + j;
  	  auto m3_row = j * k1 + i;

  	  d(m1_row,0) = m1(m1_row,0) - m3(m3_row,0);
  	}
  }

  return d;
}

int main(int argc, char* argv[])
{
#ifdef USE_SLEPC
  SlepcInitialize(&argc, &argv, NULL, NULL);
#endif

  int lam1, mu1, n;
  std::cout << "Please enter lambda1, mu1: ";
  std::cin >> lam1 >> mu1;

  std::cout << "Please enter n: ";
  std::cin >> n;

  u3::irrep u3A(lam1, mu1), u3B(n, 0), u3B_conj(0, n);

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
  u3::irrep u3C = u3Cs[index];
  int lam3 = u3C.lam(), mu3 = u3C.mu(), rhomax = rhomaxs[index];

  std::vector<u3::physicalvector> pArrA = u3::physicalBasis(u3A), 
  	pArrB = u3::physicalBasis(u3B), pArrC = u3::physicalBasis(u3C);

  double dimA = u3A.dim(), dimC = u3C.dim();

  qcx_utils::SO3_CGC cgcs;
  u3::WigU3U2_db wigU3U2db(cgcs);
  u3::WigU3SO3_db wigU3SO3db(cgcs, wigU3U2db);

  for (const auto& pA : pArrA)
  {
  	for (const auto& pB : pArrB)
  	{
  	  for (const auto& pC : pArrC)
  	  {
  	  	auto [flag, cgcU3] = wigU3SO3db.get(lam1, mu1, pA.L,
  	  										n, 0, pB.L, lam3, mu3, pC.L);

  	  	if (flag)
  	  	{
  	  	  auto [flag_conj, cgcU3_conj] = wigU3SO3db.get(lam3, mu3, pC.L,
  	  										0, n, pB.L, lam1, mu1, pA.L);

  	  	  int phi = lam1 + mu1 - pA.L - (lam3 + mu3 - pC.L) + (n + pB.L) / 2;
  	  	  double phase = (phi % 2 == 0) ? 1.0 : -1.0;
  	  	  double factor = 
  	  	  	std::sqrt((2.0 * pA.L + 1.0) * dimC / (2.0 * pC.L + 1.0) / dimA);

  	  	  std::cout << "\nl1 = " << pA.L << "; l2 = " << pB.L << "; l3 = " << pC.L 
  	  	  			<< "\nCoupling:\n" << cgcU3 
  	  	  			<< "\nConjugate coupling:\n" << cgcU3_conj 
  	  	  			<< "\nDifference with factor:\n" << diff(
  	  	  					pA.kappa, cgcU3, pC.kappa, phase * factor * cgcU3_conj)
  	  	  			<< std::endl;
  	  	}
  	  }
  	}
  }

#ifdef USE_SLEPC 
  SlepcFinalize();
#endif

  return 0;
}
