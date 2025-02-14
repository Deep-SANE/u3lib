#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <locale>
#include <string>
#include <utility>
#include <vector>

#include <decomposition.h>
#include <U3basis.h>

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

auto unique_n12n22(const std::vector<u3::gelfand>& basis)
{
  std::vector<std::pair<int, int>> ub;

  for (const auto& g : basis)
  {
    auto n12 = g.n12();
    auto n22 = g.n22();

    if (std::find(ub.cbegin(), ub.cend(), std::pair(n12, n22)) == ub.end())
      ub.emplace_back(n12, n22);
  }

  return ub;
}

int main()
{
  auto u3a = read_irrep(true, "1st");
  auto u3b = read_irrep(true, "2nd");
  auto u3c = read_irrep(false, "3rd");

  int rhomax;
  {
    auto [u3cs, rhomaxs] = u3::decomposition(u3a, u3b);

    auto iter = std::find(u3cs.cbegin(), u3cs.cend(), u3c);
    if (iter == u3cs.end())
    {
      std::cerr
          << "The 3rd irrep has not been found in the decomposition of the "
             "1st and the 2nd!"
          << std::endl;
      return EXIT_FAILURE;
    }

    auto index = std::distance(u3cs.cbegin(), iter);
    assert(u3cs[index] == u3c);

    rhomax = rhomaxs[index];
  }

  std::cout.imbue(std::locale("en_US.UTF8"));

  std::size_t NA = u3a.dim(), NB = u3b.dim(), NC = u3c.dim();

  std::cout << "Size of 1st basis: " << NA << std::endl;
  std::cout << "Size of 2nd basis: " << NB << std::endl;
  std::cout << "Size of 3rd basis: " << NC << std::endl;

  std::cout << "Total number of basis state combinations: " << NA * NB * NC
            << std::endl;

  auto t0 = std::chrono::high_resolution_clock::now();

  std::size_t reduced = 0;
  std::size_t selected = 0;

  for (int n12a = u3a.n13(); n12a >= u3a.n23(); n12a--)
    for (int n22a = u3a.n23(); n22a >= u3a.n33(); n22a--)
    {
      u3::gelfand gA(u3a.n13(), u3a.n23(), u3a.n33(), n12a, n22a, n22a);

      for (int n12b = u3b.n13(); n12b >= u3b.n23(); n12b--)
        for (int n22b = u3b.n23(); n22b >= u3b.n33(); n22b--)
        {
          u3::gelfand gB(u3b.n13(), u3b.n23(), u3b.n33(), n12b, n22b, n22b);

          for (int n12c = u3c.n13(); n12c >= u3c.n23(); n12c--)
            for (int n22c = u3c.n23(); n22c >= u3c.n33(); n22c--)
            {
              u3::gelfand gC(u3c.n13(), u3c.n23(), u3c.n33(), n12c, n22c, n22c);

              reduced++;

              bool selectionI = gA.epsilon() + gB.epsilon() == gC.epsilon();

              bool selectionII_a =
                  gC.Lambdax2() >= std::abs(gA.Lambdax2() - gB.Lambdax2());

              bool selectionII_b =
                  gC.Lambdax2() <= gA.Lambdax2() + gB.Lambdax2();

              if (selectionI && selectionII_a && selectionII_b)
                selected++;
            }
        }
    }

  auto t1 = std::chrono::high_resolution_clock::now();

  std::cout << "Total number of reduced basis state combinations: " << reduced
            << std::endl;

  std::cout << "Number of allowed reduced basis state combinations: "
            << selected << std::endl;

  std::chrono::duration<double> duration = t1 - t0;
  std::cout << "Elapsed time: " << duration.count() << " [s]" << std::endl;
}
