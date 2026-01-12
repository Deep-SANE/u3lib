#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include <u3lib.h>
#include <U3basis.h>
#include <U3SO3_to_U3U2U1.h> 

#include <qcx_utils.h>

int main(int argc, char* argv[])
{
  qcx_utils::initialize(&argc, &argv);
  u3::initialize();

#ifdef ADD_SU3LIB_PHASE_U3SO3
  int lambda, mu;
  std::cout << "Enter SU(3) quantum numbers (lambda != mu): ";
  std::cin >> lambda >> mu;

  u3::irrep ir(lambda, mu), ir_conj(mu,lambda);

  // building the gelfand basis
  std::vector<u3::gelfand> gArr = u3::gelfandBasis(ir),
                           gArr_conj = u3::gelfandBasis(ir_conj);

  // building the physical basis excluding those which have zero multiplicity
  std::vector<u3::physicalvector> pArr = u3::physicalBasis(ir);

  std::cout << "The physical basis states are: " << std::endl;
  for (int i = 0; i < pArr.size(); i++)
  {
    std::cout << "(" << i << "): |" << "(" << pArr[i].lam() << ","
              << pArr[i].mu() << "), kappa_max=" << pArr[i].kappa()
              << ", L=" << pArr[i].L() << ">" << std::endl;
  } 

  int index;
  std::cout << "Enter the index of the state that you are interested: ";
  std::cin >> index;
  const auto& p = pArr[index];

  int M;
  std::cout << "Enter projection M <= L: ";
  std::cin >> M;

  u3::physicalvector p_conj(ir_conj,p.L());
  double p_phase = ((lambda + mu - p.L() + M) % 2 == 0) ? 1.0 : -1.0;

  // Start time measurement
  auto t0 = std::chrono::high_resolution_clock::now();

  // get the transformation
  auto [canoBasis ,transformation] = u3::phys_to_cano(p, M, gArr);
  auto [canoBasis_conj, transformation_conj] = u3::phys_to_cano(p_conj, -M, gArr_conj);

  // End time measurement
  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout
      << "The Gelfand states and the expansion coefficients for kappa = 1,...,"
      << p.kappa() << " are:" << std::endl;

  for (int i = 0; i < canoBasis.size(); i++)
  {
    const auto& g = gArr[canoBasis[i]];
    auto g_conj = g.conjugate();    

    std::cout << "|" << g << "> \t" << transformation.row(i) << "\t";

    for (int j = 0; j < canoBasis_conj.size(); j++)
    {
      auto g_phase = gArr_conj[canoBasis_conj[j]].conjugate_phase();
      
      if (g_conj == gArr_conj[canoBasis_conj[j]])
      {
        auto diff = p_phase * g_phase * transformation.row(i)
                  - transformation_conj.row(j);

        if (diff.norm() > 1e-12) std::cout << "\tViolate";
        else std::cout << "\tGood";

        break;
      }
    }
    
    std::cout << std::endl;
  }

  // Calculate the duration
  std::chrono::duration<double> duration = t1 - t0;

  // Output the time in seconds
  std::cout << "Elapsed time: " << duration.count() << " seconds or "
            << duration.count() / 60.0 << " minutes" << std::endl;

#else
  std::cout << "No phase enforced" << std::endl;
#endif

  u3::finalize();
  qcx_utils::finalize();
}
