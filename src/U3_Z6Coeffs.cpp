#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string> 
#include <vector>
#include <omp.h>

#include <u3lib.h>
#include <U3decomposition.h>
#include <U3_U6Z6Coeffs.h>
#include <U3basis.h>

#include <utils.h>

int main(int argc, char* argv[])
{
  u3::initialize(&argc, &argv);

  std::cout
      << "This program compute Z6-recoupling coefficients of three U(3) irreps"
      << std::endl;

  int n13a, n23a, n33a;
  std::cout << "Please enter the I. irrep [n13a,n23a,n33a]: ";
  std::cin >> n13a >> n23a >> n33a;
  u3::irrep u3a(n13a, n23a, n33a);
  if (n33a != 0)
  {
    u3a.normalize();
    std::cout << "Your irrep is normalized to " << u3a << std::endl;
  }

  int n13b, n23b, n33b;
  std::cout << "Please enter the II. irrep [n13b,n23b,n33b]: ";
  std::cin >> n13b >> n23b >> n33b;
  u3::irrep u3b(n13b, n23b, n33b);
  if (n33b != 0)
  {
    u3b.normalize();
    std::cout << "Your irrep is normalized to " << u3b << std::endl;
  }

  int n13c, n23c, n33c;
  std::cout << "Please enter the III. irrep [n13c,n23c,n33c]: ";
  std::cin >> n13c >> n23c >> n33c;
  u3::irrep u3c(n13c, n23c, n33c);
  if (n33c != 0)
  {
    u3c.normalize();
    std::cout << "Your irrep is normalized to " << u3c << std::endl;
  }

  // couple the I. and II. irreps
  auto [u3ab, rhomaxsab] = u3::decomposition(u3a, u3b);

  std::cout << "The normalized coupled irreps of the I and II are: "
            << std::endl;
  for (int i = 0; i < u3ab.size(); i++)
  {
    std::cout << "(" << i << "): " << u3ab[i] << "\t rho = 1 ... "
              << rhomaxsab[i] << std::endl;
  }

  // couple the I. and III. irreps
  auto [u3ac, rhomaxsac] = u3::decomposition(u3a, u3c);

  std::cout << "The normalized coupled irreps of the I and III are: "
            << std::endl;
  for (int i = 0; i < u3ac.size(); i++)
  {
    std::cout << "(" << i << "): " << u3ac[i] << "\t rho = 1 ... "
              << rhomaxsac[i] << std::endl;
  }

  int i_ab, i_ac, i_ab_c;
  std::cout << "Please enter the index of the I-II coupled irrep you are "
               "interested: ";
  std::cin >> i_ab;
  std::cout << "Please enter the index of the I-III coupled irrep you are "
               "interested: ";
  std::cin >> i_ac;

  // couple the three irreps in scheme (I-II)-III
  auto [u3ab_c, rhomaxsab_c] = u3::decomposition(u3ab[i_ab], u3c);

  std::cout << "The normalized coupled irreps of the I, II and III are: "
            << std::endl;
  for (int i = 0; i < u3ab_c.size(); i++)
  {
    std::cout << "(" << i << "): " << u3ab_c[i] << "\t rho = 1 ... "
              << rhomaxsab_c[i] << std::endl;
  }
  std::cout << "Please enter the indices of the I-II-III coupled irrep you are "
               "interested: ";
  std::cin >> i_ab_c;

  // couple the three irreps in scheme (I-III)-II
  auto rhomaxac_b = u3::Rho(u3ac[i_ac], u3b, u3ab_c[i_ab_c]);

  if (rhomaxac_b > 0)
  {
    // Start time measurement
    auto t0 = std::chrono::high_resolution_clock::now();

    u3::gelfands gArrA = u3::gelfandBasis(u3a);
    u3::gelfands gArrB = u3::gelfandBasis(u3b);
    u3::gelfands gArrC = u3::gelfandBasis(u3c);
    u3::gelfands gArrAB = u3::gelfandBasis(u3ab[i_ab]);
    u3::gelfands gArrAC = u3::gelfandBasis(u3ac[i_ac]);
    u3::gelfands gArrABC = u3::gelfandBasis(u3ab_c[i_ab_c]);

    auto coeffs = u3::Z6(
        u3b, gArrB, u3a, gArrA, u3ab_c[i_ab_c], gArrABC, u3c, gArrC, u3ab[i_ab],
        gArrAB, rhomaxsab[i_ab], rhomaxsab_c[i_ab_c], u3ac[i_ac], gArrAC,
        rhomaxsac[i_ac], rhomaxac_b);

    // End time measurement
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Z([ni3b], [ni3a], [ni3],  [ni3c] "
                 ";[ni3ab],rab,rab_c,[ni3ac],rac,rac_b)"
              << std::endl;
    std::cout << "  " << u3b << "," << u3a << "," << u3ab_c[i_ab_c] << ","
              << u3c << ";" << u3ab[i_ab] << ",\t    " << u3ac[i_ac]
              << std::endl;

    std::cout << coeffs << std::endl;

    // -----------------------------------------------------------------------
    // Same coefficients via the coefficient DATABASE (u3z6db). get_U3Z6 takes
    // the six irreps (note the leading B,A ordering) and the four outer
    // multiplicities; it matches the direct Z6 result above.
    // -----------------------------------------------------------------------
    auto z6db = u3::u3z6db().get_U3Z6(u3b, u3a, u3ab_c[i_ab_c], u3c, u3ab[i_ab],
                                      rhomaxsab[i_ab], rhomaxsab_c[i_ab_c],
                                      u3ac[i_ac], rhomaxsac[i_ac], rhomaxac_b);
    std::cout << "[Database interface] u3z6db().get_U3Z6(...):\n"
              << z6db << std::endl;

    // Calculate the duration
    std::chrono::duration<double> duration = t1 - t0;

    // Output the time in seconds
    std::cout << "Elapsed time: " << duration.count() << " seconds or "
              << duration.count() / 60.0 << " minutes" << std::endl;

    // -----------------------------------------------------------------------
    // Multi-threading demonstration:
    // Each thread calls u3::initialize_thread() before database access,
    // evaluates coefficients inside, prints them out along with thread index,
    // and calls u3::finalize_thread().
    // -----------------------------------------------------------------------
    std::cout << "\n[Multi-threading demonstration]" << std::endl;
#pragma omp parallel num_threads(4)
    {
      u3::initialize_thread();
      int tid = omp_get_thread_num();

      auto blk = u3::u3z6db().get_U3Z6(u3b, u3a, u3ab_c[i_ab_c], u3c, u3ab[i_ab],
                                       rhomaxsab[i_ab], rhomaxsab_c[i_ab_c],
                                       u3ac[i_ac], rhomaxsac[i_ac], rhomaxac_b);

#pragma omp critical
      {
        std::cout << "[Thread " << tid << "] Evaluated Z6 block with size " << blk.size()
                  << ", (0,0,0,0) = " << (blk.size() > 0 ? blk(0,0,0,0) : 0.0)
                  << std::endl;
      }

      u3::finalize_thread();
    }
  }
  else
  {
    std::cout << "Recoupling does not exist" << std::endl;
  }

  u3::finalize();
}
