#ifndef UTILS_TREE_GENERATOR_H
#define UTILS_TREE_GENERATOR_H

#include <cassert>
#include <stack>
#include <vector>

namespace utils
{
  template <typename QE>
  class tree_generator
  {
    std::vector<std::stack<QE>> stacks_;

  public:
    template <typename R, typename OP, typename T>
    void generate(const R& range, OP&& op, std::vector<std::vector<T>>& result)
    {
      auto N = std::size(range);
      assert(N > 1);

      stacks_.clear();
      stacks_.resize(N - 1);

      result.clear();

      std::vector<T> row(N - 1);

      op(range[0], range[1], stacks_[0]);

      int i = 0;

      while (i >= 0)
      {
        auto& s = stacks_[i];

        while (!s.empty())
        {
          row[i] = s.top();
          s.pop();

          if (i < N - 2)
          {
            op(row[i], range[i + 2], stacks_[i + 1]);
            i++;
            break;
          }
          else
            result.push_back(row);
        }

        if (stacks_[i].empty())
          i--;
      }
    }
  };
} // namespace utils

#endif
