#include <cassert>

#include <qcx_print.h>

#include <U3irrep.h>

using qcx_utils::print;

int main()
{
  su3::irrep su3A;
  su3::irrep su3B(2, 1);

  print("{}\n", su3A);
  print("{}\n", su3B);
}
