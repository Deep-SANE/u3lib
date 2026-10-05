#pragma once
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <u3lib.h>

namespace test
{
  inline void require(bool ok, const std::string& what)
  {
    if (!ok)
      throw std::runtime_error(what);
  }
  inline void near(double actual, double expected, const std::string& what,
                   double atol = 1e-12, double rtol = 1e-10)
  {
    require(std::isfinite(actual) && std::isfinite(expected) &&
                std::abs(actual - expected) <= atol + rtol * std::abs(expected),
            what);
  }
  inline double sign(int n) { return n % 2 ? -1.0 : 1.0; }
  inline std::string irname(const u3::irrep& ir)
  {
    return "(" + std::to_string(ir.lam()) + "," + std::to_string(ir.mu()) + ")";
  }
  template <class... T>
  std::string key(T... args)
  {
    std::ostringstream out;
    ((out << args << ':'), ...);
    return out.str();
  }
  // Generation is a separate executable, never enabled by the test runner.
  class Reference
  {
    std::map<std::string, double> expected;
    std::ofstream output;
    std::string filename;

  public:
    Reference(const std::string& name,
              const std::string& directory = "reference", bool generate = false)
        : filename(directory + "/" + name + ".dat")
    {
      if (generate)
      {
        output.open(filename);
        require(bool(output), "Cannot write " + filename);
        output << "# Signed RREF-canonicalized reference values; labels and "
                  "provenance: reference/README.md\n"
               << std::setprecision(17);
      }
      else
      {
        std::ifstream in(filename);
        require(bool(in), "Cannot read " + filename);
        std::string line;
        while (std::getline(in, line))
        {
          if (line.empty() || line[0] == '#')
            continue;
          std::istringstream row(line);
          std::string label, extra;
          double value;
          require(bool(row >> label >> value) && !(row >> extra) &&
                      std::isfinite(value),
                  "Malformed reference in " + filename);
          require(expected.emplace(label, value).second,
                  "Duplicate reference " + label);
        }
        require(!expected.empty(), "Empty reference " + filename);
      }
    }
    void check(const std::string& label, double value)
    {
      require(std::isfinite(value), "Nonfinite result: " + label);
      if (output.is_open())
      {
        output << label << ' ' << value << '\n';
        return;
      }
      auto it = expected.find(label);
      require(it != expected.end(), "Missing documented value: " + label);
      near(value, it->second, "Documented-value mismatch: " + label);
      expected.erase(it);
    }
    void finish()
    {
      require(expected.empty(), "Unused documented values in " + filename);
      if (output.is_open())
      {
        output.flush();
        require(bool(output), "Write failed: " + filename);
      }
    }
  };
  struct Suite
  {
    int passed = 0, failed = 0;
    template <class F>
    void run(const std::string& label, F work)
    {
      std::cout << "[TEST] " << label << std::endl;
      try
      {
        work();
        ++passed;
        std::cout << "[PASS] " << label << std::endl;
      }
      catch (const std::exception& e)
      {
        ++failed;
        std::cout << "[FAIL] " << label << ": " << e.what() << std::endl;
      }
    }
    int result() const
    {
      std::cout << "Summary: " << passed << " passed, " << failed
                << " failed\n";
      return failed ? 1 : 0;
    }
  };
  inline Eigen::MatrixXd full(const std::vector<int>& ids,
                              const Eigen::MatrixXd& v, int dim)
  {
    Eigen::MatrixXd f = Eigen::MatrixXd::Zero(dim, v.cols());
    for (int i = 0; i < int(ids.size()); ++i)
      f.row(ids[i]) = v.row(i);
    return f;
  }
  inline void orthogonal(const Eigen::MatrixXd& v, const std::string& label)
  {
    near((v.transpose() * v - Eigen::MatrixXd::Identity(v.cols(), v.cols()))
             .norm(),
         0, label, 1e-9, 0);
  }
  struct Coupling
  {
    u3::irrep a, b, c;
  };
  inline std::vector<Coupling> couplings()
  {
    return {{{0, 0}, {5, 3}, {5, 3}}, {{1, 0}, {1, 0}, {2, 0}},
            {{1, 0}, {1, 0}, {0, 1}}, {{1, 0}, {0, 1}, {0, 0}},
            {{1, 0}, {0, 1}, {1, 1}}, {{1, 1}, {1, 1}, {1, 1}},
            {{2, 2}, {1, 1}, {2, 2}}, {{3, 2}, {2, 3}, {2, 2}},
            {{4, 4}, {1, 1}, {4, 4}}, {{5, 3}, {1, 1}, {5, 3}},
            {{5, 5}, {1, 1}, {5, 5}}, {{5, 0}, {0, 5}, {0, 0}}};
  }
  inline std::string label(const Coupling& x)
  {
    return irname(x.a) + "x" + irname(x.b) + "->" + irname(x.c);
  }
} // namespace test
