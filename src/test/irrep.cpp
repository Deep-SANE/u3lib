#include <cassert>

#include <qcx_print.h>

#include <U3irrep.h>

using qcx_utils::print;

void print_irrep(const auto& i)
{
  print("{}: dim = {}\n", i, i.dim());
}

int main()
{
  su3::irrep su3A;
  su3::irrep su3B(2, 1);
  assert(su3B != su3A);

  print_irrep(su3A);
  print_irrep(su3B);

  su3::irrep su3C = su3B;
  assert(su3B == su3C);

  u3::irrep u3A = su3A;
  u3::irrep u3B = su3B;
  assert(u3A != u3B);

  print_irrep(u3A);
  print_irrep(u3B);

  u3::irrep u3C = u3B;
  assert(u3B == u3C);

  u3::irrep u3D(5,3,2);
  print_irrep(u3D);

  su3::irrep su3D = u3D;
  print_irrep(su3D);

  u3D.normalize();
  print_irrep(u3D);

  su3D = u3D;
  print_irrep(su3D);

  {
    su3::irrep su3(2, 2);
    u3::irrep u3(su3, 6);
    print_irrep(u3);

    assert(u3.n13() == 4);
    assert(u3.n23() == 2);
    assert(u3.n33() == 0);
  }

  {
    su3::irrep su3(0, 0);
    u3::irrep u3(su3, 6);
    print_irrep(u3);

    assert(u3.n13() == 2);
    assert(u3.n23() == 2);
    assert(u3.n33() == 2);
  }
}
