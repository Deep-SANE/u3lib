#include <chrono>
#include <fstream> 
#include <iomanip>
#include <iostream>
#include <sstream> 
#include <string>
#include <vector> 

#include <u3lib.h>
#include <U3decomposition.h>
#include <U3_U9Coeffs.h> 
#include <U3basis.h>

#include <qcx_utils.h>

int main(int argc, char* argv[])
{
  qcx_utils::initialize(&argc, &argv);
  u3::initialize();

  std::cout
      << "This program compute U9-recoupling coefficients of four U(3) irreps"
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

  int n13d, n23d, n33d;
  std::cout << "Please enter the IV. irrep [n13d,n23d,n33d]: ";
  std::cin >> n13d >> n23d >> n33d;
  u3::irrep u3d(n13d, n23d, n33d);
  if (n33d != 0)
  {
    u3d.normalize();
    std::cout << "Your irrep is normalized to " << u3d << std::endl;
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
    std::cout << "(" << i << "): " << u3ac[i]<< "\t rho = 1 ... "
              << rhomaxsac[i] << std::endl;
  }

  // couple the II. and IV. irreps
  auto [u3bd, rhomaxsbd] = u3::decomposition(u3b, u3d);

  std::cout << "The normalized coupled irreps of the II and IV are: "
            << std::endl;
  for (int i = 0; i < u3bd.size(); i++)
  {
    std::cout << "(" << i << "): " << u3bd[i] << "\t rho = 1 ... "
              << rhomaxsbd[i] << std::endl;
  }

  // couple the III. and IV. irreps
  auto [u3cd, rhomaxscd] = u3::decomposition(u3c, u3d);

  std::cout << "The normalized coupled irreps of the III and IV are: "
            << std::endl;
  for (int i = 0; i < u3cd.size(); i++)
  {
    std::cout << "(" << i << "): " << u3cd[i] << "\t rho = 1 ... "
              << rhomaxscd[i] << std::endl;
  }

  int i_ab, i_ac, i_bd, i_cd, i_ab_cd, i_ac_bd;
  std::cout << "Please enter the index of the I-II coupled irrep you are "
               "interested: ";
  std::cin >> i_ab;
  std::cout << "Please enter the index of the I-III coupled irrep you are "
               "interested: ";
  std::cin >> i_ac;
  std::cout << "Please enter the index of the II-IV coupled irrep you are "
               "interested: ";
  std::cin >> i_bd;
  std::cout << "Please enter the index of the III-IV coupled irrep you are "
               "interested: ";
  std::cin >> i_cd;

  // couple the 4 irreps in the scheme (I-II)-(III-IV)
  auto [u3ab_cd, rhomaxsab_cd] = u3::decomposition(u3ab[i_ab], u3cd[i_cd]);

  std::cout << "The coupled irreps of the I, II, III and IV are: " << std::endl;
  for (int i = 0; i < u3ab_cd.size(); i++)
  {
    std::cout << "(" << i << "): " << u3ab_cd[i] << "\t rho = 1 ... " 
              << rhomaxsab_cd[i] << std::endl;
  }
  std::cout << "Please enter the index of the I-II-III-IV coupled irrep you "
               "are interested: ";
  std::cin >> i_ab_cd;

  // couple the 4 irreps in scheme (I-III)-(II-IV)
  auto rhomaxac_bd = u3::Rho(u3ac[i_ac], u3bd[i_bd], u3ab_cd[i_ab_cd]);

  if (rhomaxac_bd > 0)
  {
    // Start time measurement
    auto t0 = std::chrono::high_resolution_clock::now();

    u3::gelfands gArrA = u3::gelfandBasis(u3a);
    u3::gelfands gArrB = u3::gelfandBasis(u3b);
    u3::gelfands gArrC = u3::gelfandBasis(u3c);
    u3::gelfands gArrD = u3::gelfandBasis(u3d);
    u3::gelfands gArrAB = u3::gelfandBasis(u3ab[i_ab]);
    u3::gelfands gArrCD = u3::gelfandBasis(u3cd[i_cd]);
    u3::gelfands gArrAC = u3::gelfandBasis(u3ac[i_ac]);
    u3::gelfands gArrBD = u3::gelfandBasis(u3bd[i_bd]);
    u3::gelfands gArrABCD = u3::gelfandBasis(u3ab_cd[i_ab_cd]);

    auto coeffs =
      u3::U9(u3a, gArrA, u3b, gArrB, u3c, gArrC, u3d, gArrD, u3ab[i_ab],
             gArrAB, u3cd[i_cd], gArrCD, u3ac[i_ac], gArrAC, u3bd[i_bd],
             gArrBD, u3ab_cd[i_ab_cd], gArrABCD, rhomaxsab[i_ab],
             rhomaxscd[i_cd], rhomaxsab_cd[i_ab_cd], rhomaxsac[i_ac],
             rhomaxsbd[i_bd], rhomaxac_bd);

    // End time measurement
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << u3a << "\t" << u3b << "\t" << u3ab[i_ab] << "\trhoAB"
              << std::endl;
    std::cout << u3c << "\t" << u3d << "\t" << u3cd[i_cd] << "\trhoCD"
              << std::endl;
    std::cout << u3ac[i_ac] << "\t" << u3bd[i_bd] << "\t" << u3ab_cd[i_ab_cd]
              << "\trhoACBD" << std::endl;
    std::cout << "rhoAC\trhoBD\trhoABCD\t" << std::endl;

    std::cout << coeffs << std::endl;
    
    // Calculate the duration
    std::chrono::duration<double> duration = t1 - t0;

    // Output the time in seconds
    std::cout << "Elapsed time: " << duration.count() << " seconds or "
              << duration.count() / 60.0 << " minutes" << std::endl;
  }
  else
  {
    std::cout << "Recoupling does not exist" << std::endl;
  }

  u3::finalize();
  qcx_utils::finalize();
}
