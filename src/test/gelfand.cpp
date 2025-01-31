#ifdef NDEBUG
  #error "Test compiled with disabled assertions; do not define NDEBUG!"
#endif

#include <cassert>
#include <iostream>

#include <U3basis.h>

int main()
{
  // default constructor:
  u3::gelfand_t g1;

  assert(g1.n13() == 0);
  assert(g1.n23() == 0);
  assert(g1.n33() == 0);
  assert(g1.n12() == 0);
  assert(g1.n22() == 0);
  assert(g1.n11() == 0);

  // converting constructor:
  u3::gelfand_t g2(0, 0, 0, 0, 0, 0);

  assert(g2.n13() == 0);
  assert(g2.n23() == 0);
  assert(g2.n33() == 0);
  assert(g2.n12() == 0);
  assert(g2.n22() == 0);
  assert(g2.n11() == 0);

  // comparison operators (equal gps):
  assert(g2 == g1);
  assert(!(g2 != g1));

  u3::gelfand_t g3(8, 4, 0, 6, 2, 3);

  assert(g3.n13() == 8);
  assert(g3.n23() == 4);
  assert(g3.n33() == 0);
  assert(g3.n12() == 6);
  assert(g3.n22() == 2);
  assert(g3.n11() == 3);

  // comparision operators (different gps):
  assert(g3 != g1);
  assert(!(g3 == g1));

  // normalization:
  u3::gelfand_t g4(9, 5, 1, 7, 3, 4);
  g4.normalize();

  assert(g4.n13() == 8);
  assert(g4.n23() == 4);
  assert(g4.n33() == 0);
  assert(g4.n12() == 6);
  assert(g4.n22() == 2);
  assert(g4.n11() == 3);

  assert(g4 == g3);

  // setters:
  u3::gelfand_t g5;

  g5.n13(8);
  g5.n23(4);
  g5.n33(0);
  g5.n12(6);
  g5.n22(2);
  g5.n11(3);

  assert(g5.n13() == 8);
  assert(g5.n23() == 4);
  assert(g5.n33() == 0);
  assert(g5.n12() == 6);
  assert(g5.n22() == 2);
  assert(g5.n11() == 3);

  assert(g5 == g4);

  // copy constructor:
  u3::gelfand_t g6 = g5;

  assert(g6.n13() == 8);
  assert(g6.n23() == 4);
  assert(g6.n33() == 0);
  assert(g6.n12() == 6);
  assert(g6.n22() == 2);
  assert(g6.n11() == 3);

  assert(g6 == g5);

  // copy assignment:
  u3::gelfand_t g7;
  g7 = g6;

  assert(g7.n13() == 8);
  assert(g7.n23() == 4);
  assert(g7.n33() == 0);
  assert(g7.n12() == 6);
  assert(g7.n22() == 2);
  assert(g7.n11() == 3);

  assert(g7 == g6);

  std::cout << "Seemingly, all tests have passed." << std::endl;
}
