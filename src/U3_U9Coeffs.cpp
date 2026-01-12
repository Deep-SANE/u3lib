#include <chrono>
#include <fstream> 
#include <iomanip>
#include <iostream>
#include <sstream> 
#include <string>
#include <vector> 

#include <U3decomposition.h>
#include <U3_U9Coeffs.h> 
#include <U3basis.h>

#include <qcx_utils.h>

int main(int argc, char* argv[])
{
  qcx_utils::initialize(&argc, &argv);

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

    std::vector<u3::gelfand> gArrA = u3::gelfandBasis(u3a);
    std::vector<u3::gelfand> gArrB = u3::gelfandBasis(u3b);
    std::vector<u3::gelfand> gArrC = u3::gelfandBasis(u3c);
    std::vector<u3::gelfand> gArrD = u3::gelfandBasis(u3d);
    std::vector<u3::gelfand> gArrAB = u3::gelfandBasis(u3ab[i_ab]);
    std::vector<u3::gelfand> gArrCD = u3::gelfandBasis(u3cd[i_cd]);
    std::vector<u3::gelfand> gArrAC = u3::gelfandBasis(u3ac[i_ac]);
    std::vector<u3::gelfand> gArrBD = u3::gelfandBasis(u3bd[i_bd]);
    std::vector<u3::gelfand> gArrABCD = u3::gelfandBasis(u3ab_cd[i_ab_cd]);

    auto& cgcs = qcx_utils::cgcs();
    auto& racahs = qcx_utils::racahs();
    u3::WigU3U2_db wigU3U2db;

    std::vector<double> coeffs =
      u3::U9(u3a, gArrA, u3b, gArrB, u3c, gArrC, u3d, gArrD, u3ab[i_ab],
             gArrAB, u3cd[i_cd], gArrCD, u3ac[i_ac], gArrAC, u3bd[i_bd],
             gArrBD, u3ab_cd[i_ab_cd], gArrABCD, rhomaxsab[i_ab],
             rhomaxscd[i_cd], rhomaxsab_cd[i_ab_cd], rhomaxsac[i_ac],
             rhomaxsbd[i_bd], rhomaxac_bd, racahs, wigU3U2db);

    // End time measurement
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << u3a << "\t" << u3b << "\t" << u3ab[i_ab] << "\trhoab"
              << std::endl;
    std::cout << u3c << "\t" << u3d << "\t" << u3cd[i_cd] << "\trhocd"
              << std::endl;
    std::cout << u3ac[i_ac] << "\t" << u3bd[i_bd] << "\t" << u3ab_cd[i_ab_cd]
              << "\trhoac_bd" << std::endl;
    std::cout << "rhoac\trhobd\trhoab_cd\t" << std::endl;

    std::cout << std::endl;

    std::cout << "rhoac\trhobd\trhoab_cd\trhoab\trhocd\trhoac_bd" << std::endl;

    for (int rac = 0; rac < rhomaxsac[i_ac]; rac++)
      for (int rbd = 0; rbd < rhomaxsbd[i_bd]; rbd++)
        for (int rab_cd = 0; rab_cd < rhomaxsab_cd[i_ab_cd]; rab_cd++)
          for (int rab = 0; rab < rhomaxsab[i_ab]; rab++)
            for (int rcd = 0; rcd < rhomaxscd[i_cd]; rcd++)
              for (int rac_bd = 0; rac_bd < rhomaxac_bd; rac_bd++)
                std::cout << "  " << rac + 1 << "\t  " << rbd + 1 << "\t  "
                          << rab_cd + 1 << "\t\t  " << rab + 1 << "\t  "
                          << rcd + 1 << "\t  " << rac_bd + 1 << "\t\t  ->  "
                          << coeffs[((((rac * rhomaxsbd[i_bd] + rbd) *
                                        rhomaxsab_cd[i_ab_cd] + rab_cd) *
                                        rhomaxsab[i_ab] + rab) *
                                        rhomaxscd[i_cd] + rcd) *
                                        rhomaxac_bd + rac_bd]
                          << std::endl;

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

  qcx_utils::finalize();
}
