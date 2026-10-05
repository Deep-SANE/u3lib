#include "conjugation_cases.h"
int main(int argc, char** argv)
{
  u3::initialize(&argc, &argv);
  test::Suite suite;
  try
  {
    test::Reference ref("conjugation");
    test::conjugation(suite, ref);
    suite.run("conjugation: all documented entries exercised",
              [&] { ref.finish(); });
  }
  catch (const std::exception& e)
  {
    ++suite.failed;
    std::cout << "[FAIL] " << e.what() << '\n';
  }
  u3::finalize();
  return suite.result();
}
