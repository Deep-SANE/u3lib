#ifndef U3LIB_SO3_RACAH_H
#define U3LIB_SO3_RACAH_H

#include <cmath>

extern "C"
{
#include <wigxjpf.h>
}

namespace u3
{
  class SO3_Racah
  {
  public:
    SO3_Racah()
    {
      wig_table_init(2 * 100, 9);
      wig_temp_init(2 * 100);
    }

    ~SO3_Racah()
    {
      wig_temp_free();
      wig_table_free();
    }

    double compute(int ll1, int ll2, int ll, int ll3, int ll12, int ll23)
    {
      double phase = ((ll1 + ll2 + ll3 + ll) % 4 == 0) ? 1 : -1;
      double factor = std::sqrt((static_cast<double>(ll12) + 1.0) * (static_cast<double>(ll23) + 1.0));
      double val6j = phase * factor * wig6jj(ll1, ll2, ll12, ll3, ll, ll23);

      return val6j;
    }
  };
} // namespace u3

#endif
