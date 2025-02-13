#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <locale>
#include <string>
#include <utility>
#include <vector>

#include <decomposition.h>
#include <U3basis.h>
#include <U3U2U1xU3U2U1_2.h>

auto read_irrep(bool normalize, const std::string& nth = "")
{
  int n13, n23, n33;

  std::cout << "Please enter the " << nth << " irrep [n13, n23, n33]: ";
  std::cin >> n13 >> n23 >> n33;

  u3::irrep irrep(n13, n23, n33);

  if (normalize && n33)
  {
    irrep.normalize();
    std::cout << "Your irrep has been normalized to " << irrep << std::endl;
  }

  std::cout << std::endl;

  return irrep;
}

auto unique_SU2_basis(const std::vector<u3::gelfand>& basis)
{
  std::vector<std::pair<int, int>> ub;

  for (const auto& g : basis)
  {
    auto n12 = g.n12() - g.n33();
    auto n22 = g.n22() - g.n33();

    if (std::find(ub.cbegin(), ub.cend(), std::pair(n12, n22)) == ub.end())
      ub.emplace_back(n12, n22);
  }

  return ub;
}

int main()
{
  auto u3a = read_irrep(true, "1st");
  auto u3b = read_irrep(true, "2nd");

  auto [u3cs, rhomaxs] = u3::decomposition(u3a, u3b);

  auto u3c = read_irrep(false, "3rd");

  auto iter = std::find(u3cs.cbegin(), u3cs.cend(), u3c);
  if (iter == u3cs.end())
  {
    std::cerr << "The 3rd irrep has not been found in the decomposition of the "
                 "1st and the 2nd!"
              << std::endl;
    return EXIT_FAILURE;
  }

  auto index = std::distance(u3cs.cbegin(), iter);
  assert(u3cs[index] == u3c);

  auto rhomax = rhomaxs[index];

  auto n33c = u3c.n33();

  auto gArrA = u3::gelfandBasis(u3a);
  auto gArrB = u3::gelfandBasis(u3b);
  auto gArrC = u3::gelfandBasis(u3c);

  std::cout.imbue(std::locale("en_US.UTF8"));

  std::cout << "Sizes of Gelfand bases: " << gArrA.size() << ", "
            << gArrB.size() << "; " << gArrC.size() << std::endl;
  std::cout << "Total number of basis state combinations: "
            << gArrA.size() * gArrB.size() * gArrC.size() << std::endl;

  auto reduced_gArrA = unique_SU2_basis(gArrA);
  auto reduced_gArrB = unique_SU2_basis(gArrB);
  auto reduced_gArrC = unique_SU2_basis(gArrC);

  std::cout << "Sizes of reduced Gelfand bases: " << reduced_gArrA.size()
            << ", " << reduced_gArrB.size() << "; " << reduced_gArrC.size()
            << std::endl;
  std::cout << "Total number of reduced basis state combinations: "
            << reduced_gArrA.size() * reduced_gArrB.size() *
                   reduced_gArrC.size()
            << std::endl;
}
