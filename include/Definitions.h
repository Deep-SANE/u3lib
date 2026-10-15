#ifndef U3LIB_GLOBAL_DEFINITIONS_H
#define U3LIB_GLOBAL_DEFINITIONS_H

#include <iostream>

#include <ND.h>

namespace u3
{

//////////// COMPILER DEFINITIONS //////////////

#define U3LIB_ENABLE_WIGU3U2_CACHE
// #define TIME_WIGU3U2

#define U3LIB_ENABLE_WIGU3SO3_CACHE
// #define TIME_WIGU3SO3

#define U3LIB_ENABLE_U3U6_CACHE
#define U3LIB_ENABLE_U3Z6_CACHE
// #define TIME_U3U6Z6

#define U3LIB_ENABLE_U3U9_CACHE
// #define TIME_U3U9

//////////// INTEGER DEFINITIONS //////////////

  using uqn_t = uint8_t;  // unsigned quantum numbers: 0 to 255
  using qn_t = int8_t;    // signed quantum numbers: -128 to 127
  using mult_t = uint8_t; // multiplicities: 0 to 255

//////////// COUPLING & RECOUPLING COEFFICIENTS //////////////

/* General rule
1. operator[] -> accessing using index
2. operator() -> accessing using multiplicities
3. access for non-const object: allowing modification
4. access for const object: forbidding modification
5. data() -> exposing internal data
*/
  using coeffs_t = std::vector<double>;

  template <typename T>
  class U3U2
  {
  	T coeffs_;

  public:
  	U3U2() = default;
  	U3U2(size_t rhomax) : coeffs_(rhomax, 0.0) {}
  	U3U2(T&& coeffs) : coeffs_(std::move(coeffs)) {}

  	size_t size() const { return coeffs_.size(); }
  	bool empty() const { return coeffs_.empty(); }
  	const auto& data() const { return coeffs_; }

  	size_t rhomax() const { return coeffs_.size(); }

  	double& operator[](size_t i)
  	{
  	  assert(i < coeffs_.size());
      return coeffs_[i];
  	}

  	double operator[](size_t i) const
  	{
  	  assert(i < coeffs_.size());
      return coeffs_[i];
  	}

  	double& operator()(size_t rho)
  	{
  	  return operator[](rho);
  	}

  	double operator()(size_t rho) const
  	{
  	  return operator[](rho);
  	}

    U3U2& operator*=(double s) 
    {
      for (auto& x : coeffs_) x *= s;
      return *this;
    }

    U3U2 operator*(double s) const 
    {
      U3U2 result = *this;
      result *= s;
      return result;
    }

    friend U3U2 operator*(double s, const U3U2& t) 
    {
      return t * s;
    }

    U3U2& operator/=(double s) 
    {
      assert(s != 0.0);
      for (auto& x : coeffs_) x /= s;
      return *this;
    }

    U3U2 operator/(double s) const 
    {
      assert(s != 0.0);
      U3U2 result = *this;
      result /= s;
      return result;
    }

  	friend std::ostream& operator<<(std::ostream& s, const U3U2& x)
  	{
  	  s << "rho = 0..." << (static_cast<int>(x.rhomax()) - 1) << " -> ";
  	  for (const auto& c : x.data()) s << c << "\t";
  	  return s;
  	}

    void print(std::ostream& s = std::cout) const
    {
      s << *this << "\n";
    }
  };

  template <typename T>
  class U3SO3
  {
    T coeffs_;
    utils::NDIndex4 ndi_;

  public:
    U3SO3()
    {
      ndi_.add_dimensions(0, 0, 0, 0);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    U3SO3(size_t k1max, size_t k2max, size_t k3max, size_t rhomax)
      : coeffs_(k1max * k2max * k3max * rhomax, 0.0)
    {
      ndi_.add_dimensions(k1max, k2max, k3max, rhomax); 
    }

    U3SO3(T&& coeffs,
    	  size_t k1max, size_t k2max, size_t k3max, size_t rhomax)
      : coeffs_(std::move(coeffs))
    {
      ndi_.add_dimensions(k1max, k2max, k3max, rhomax);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    size_t size() const { return coeffs_.size(); }
  	bool empty() const { return coeffs_.empty(); }
  	const auto& data() const { return coeffs_; }

    size_t kappaAmax() const { return ndi_.dimension(0); }
    size_t kappaBmax() const { return ndi_.dimension(1); }
    size_t kappaCmax() const { return ndi_.dimension(2); }
    size_t rhomax() const { return ndi_.dimension(3); }

    double& operator[](size_t i)
    {
      assert(i < coeffs_.size());
      return coeffs_[i];
    }

    double operator[](size_t i) const
    {
      assert(i < coeffs_.size());
      return coeffs_[i];
    }

    double& operator()(size_t k1, size_t k2, size_t k3, size_t rho)
    {
      std::array multi_index = {k1, k2, k3, rho};
      size_t i = ndi_.to_linear(multi_index);
      return operator[](i);
    }

    double operator()(size_t k1, size_t k2, size_t k3, size_t rho) const
    {
      std::array multi_index = {k1, k2, k3, rho};
      size_t i = ndi_.to_linear(multi_index);
      return operator[](i);
    }

    U3SO3& operator*=(double s) 
    {
      for (auto& x : coeffs_) x *= s;
      return *this;
    }

    U3SO3 operator*(double s) const 
    {
      U3SO3 result = *this;
      result *= s;
      return result;
    }

    friend U3SO3 operator*(double s, const U3SO3& t) 
    {
      return t * s;
    }

    U3SO3& operator/=(double s) 
    {
      assert(s != 0.0);
      for (auto& x : coeffs_) x /= s;
      return *this;
    }

    U3SO3 operator/(double s) const 
    {
      assert(s != 0.0);
      U3SO3 result = *this;
      result /= s;
      return result;
    }

    friend std::ostream& operator<<(std::ostream& s, const U3SO3& x)
  	{
  	  s << "kappaA\tkappaB\tkappaC\trho\n";
  	  std::vector<size_t> multi(4);
  	  for (size_t i = 0; i < x.size(); i++)
  	  {
  	    x.ndi_.to_multi(i, multi);
  	    s << multi[0] << "\t" << multi[1] << "\t"
  	      << multi[2] << "\t" << multi[3] << "\t -> "
  	      << x[i] << "\n";
  	  }
	    return s;
  	}

    void print(std::ostream& s = std::cout) const
    {
      s << *this;
    }
  };

  template <typename T>
  class U3U6
  {
    T coeffs_;
    utils::NDIndex4 ndi_;

  public:
    U3U6()
    {
      ndi_.add_dimensions(0, 0, 0, 0);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    U3U6(size_t rhoABmax, size_t rhoAB_Cmax, 
         size_t rhoBCmax, size_t rhoA_BCmax)
      : coeffs_(rhoABmax * rhoAB_Cmax * rhoBCmax * rhoA_BCmax, 0.0)
    {
      ndi_.add_dimensions(rhoABmax, rhoAB_Cmax, rhoBCmax, rhoA_BCmax);
    }

    U3U6(T&& coeffs, size_t rhoABmax, size_t rhoAB_Cmax, 
                     		size_t rhoBCmax, size_t rhoA_BCmax)
      : coeffs_(std::move(coeffs))
    {
      ndi_.add_dimensions(rhoABmax, rhoAB_Cmax, rhoBCmax, rhoA_BCmax);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    size_t size() const { return coeffs_.size(); }
  	bool empty() const { return coeffs_.empty(); }
  	const auto& data() const { return coeffs_; }

    size_t rhoABmax() const { return ndi_.dimension(0); }
    size_t rhoAB_Cmax() const { return ndi_.dimension(1); }
    size_t rhoBCmax() const { return ndi_.dimension(2); }
    size_t rhoA_BCmax() const { return ndi_.dimension(3); }

    double& operator[](size_t i)
    {
      assert(i < coeffs_.size());
      return coeffs_[i];
    }

    double operator[](size_t i) const
    {
      assert(i < coeffs_.size());
      return coeffs_[i];
    }

    double& operator()(size_t rhoAB, size_t rhoAB_C,
    				           size_t rhoBC, size_t rhoA_BC)
    {
      std::array multi_index = {rhoAB, rhoAB_C, rhoBC, rhoA_BC};
      size_t i = ndi_.to_linear(multi_index);
      return operator[](i);
    }

    double operator()(size_t rhoAB, size_t rhoAB_C,
    				          size_t rhoBC, size_t rhoA_BC) const
    {
      std::array multi_index = {rhoAB, rhoAB_C, rhoBC, rhoA_BC};
      size_t i = ndi_.to_linear(multi_index);
      return operator[](i);
    }

    friend std::ostream& operator<<(std::ostream& s, const U3U6& x)
  	{
  	  s << "rhoAB\trhoAB_C\trhoBC\trhoA_BC\n";
  	  std::vector<size_t> multi(4);
  	  for (size_t i = 0; i < x.size(); i++)
  	  {
  	    x.ndi_.to_multi(i, multi);
  	    s << multi[0] << "\t" << multi[1] << "\t"
  	      << multi[2] << "\t" << multi[3] << "\t -> "
  	      << x[i] << "\n";
  	  }
	    return s;
  	}

    void print(std::ostream& s = std::cout) const
    {
      s << *this;
    }
  };

  template <typename T>
  class U3Z6
  {
    T coeffs_;
    utils::NDIndex4 ndi_;

  public:
    U3Z6()
    {
      ndi_.add_dimensions(0, 0, 0, 0);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    U3Z6(size_t rhoABmax, size_t rhoAB_Cmax, 
         size_t rhoACmax, size_t rhoAC_Bmax)
      : coeffs_(rhoABmax * rhoAB_Cmax * rhoACmax * rhoAC_Bmax, 0.0)
    {
      ndi_.add_dimensions(rhoABmax, rhoAB_Cmax, rhoACmax, rhoAC_Bmax);
    }

    U3Z6(T&& coeffs, size_t rhoABmax, size_t rhoAB_Cmax, 
                 		    size_t rhoACmax, size_t rhoAC_Bmax)
        : coeffs_(std::move(coeffs))
    {
      ndi_.add_dimensions(rhoABmax, rhoAB_Cmax, rhoACmax, rhoAC_Bmax);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    size_t size() const { return coeffs_.size(); }
  	bool empty() const { return coeffs_.empty(); }
  	const auto& data() const { return coeffs_; }

    size_t rhoABmax() const { return ndi_.dimension(0); }
    size_t rhoAB_Cmax() const { return ndi_.dimension(1); }
    size_t rhoACmax() const { return ndi_.dimension(2); }
    size_t rhoAC_Bmax() const { return ndi_.dimension(3); }

    double& operator[](size_t i)
    {
      assert(i < coeffs_.size());
      return coeffs_[i];
    }

    double operator[](size_t i) const
    {
      assert(i < coeffs_.size());
      return coeffs_[i];
    }

    double& operator()(size_t rhoAB, size_t rhoAB_C,
    				           size_t rhoAC, size_t rhoAC_B)
    {
      std::array multi_index = {rhoAB, rhoAB_C, rhoAC, rhoAC_B};
      size_t i = ndi_.to_linear(multi_index);
      return operator[](i);
    }

    double operator()(size_t rhoAB, size_t rhoAB_C,
    				          size_t rhoAC, size_t rhoAC_B) const
    {
      std::array multi_index = {rhoAB, rhoAB_C, rhoAC, rhoAC_B};
      size_t i = ndi_.to_linear(multi_index);
      return operator[](i);
    }

    friend std::ostream& operator<<(std::ostream& s, const U3Z6& x)
  	{
  	  s << "rhoAB\trhoAB_C\trhoAC\trhoAC_B\n";
  	  std::vector<size_t> multi(4);
  	  for (size_t i = 0; i < x.size(); i++)
  	  {
  	    x.ndi_.to_multi(i, multi);
  	    s << multi[0] << "\t" << multi[1] << "\t"
  	      << multi[2] << "\t" << multi[3] << "\t -> "
  	      << x[i] << "\n";
  	  }
	    return s;
  	}

    void print(std::ostream& s = std::cout) const
    {
      s << *this;
    }
  };

  template <typename T>
  class U3U9
  {
    T coeffs_;
    utils::NDIndex6 ndi_;

  public:
    U3U9()
    {
      ndi_.add_dimensions(0, 0, 0, 0, 0, 0);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    U3U9(size_t rhoACmax, size_t rhoBDmax, size_t rhoABCDmax,
         size_t rhoABmax, size_t rhoCDmax, size_t rhoACBDmax)
        : coeffs_(rhoACmax * rhoBDmax * rhoABCDmax * 
        		  rhoABmax * rhoCDmax * rhoACBDmax, 0.0)
    {
      ndi_.add_dimensions(rhoACmax, rhoBDmax, rhoABCDmax, 
      					  rhoABmax, rhoCDmax, rhoACBDmax);
    }

    U3U9(T&& coeffs,
    	 size_t rhoACmax, size_t rhoBDmax, size_t rhoABCDmax,
         size_t rhoABmax, size_t rhoCDmax, size_t rhoACBDmax)
        : coeffs_(std::move(coeffs))
    {
      ndi_.add_dimensions(rhoACmax, rhoBDmax, rhoABCDmax, rhoABmax, rhoCDmax,
                          rhoACBDmax);
      assert(coeffs_.size() == ndi_.cardinality());
    }

    size_t size() const { return coeffs_.size(); }
  	bool empty() const { return coeffs_.empty(); }
  	const auto& data() const { return coeffs_; }

    size_t rhoACmax() const { return ndi_.dimension(0); }
    size_t rhoBDmax() const { return ndi_.dimension(1); }
    size_t rhoABCDmax() const { return ndi_.dimension(2); }
    size_t rhoABmax() const { return ndi_.dimension(3); }
    size_t rhoCDmax() const { return ndi_.dimension(4); }
    size_t rhoACBDmax() const { return ndi_.dimension(5); }

    double& operator[](size_t i)
    {
      assert(i < coeffs_.size());
      return coeffs_[i];
    }

    double operator[](size_t i) const
    {
      assert(i < coeffs_.size());
      return coeffs_[i];
    }

    double& operator()(size_t rhoAC, size_t rhoBD, size_t rhoABCD, 
    				           size_t rhoAB, size_t rhoCD, size_t rhoACBD)
    {
      std::array multi_index = {rhoAC, rhoBD, rhoABCD, rhoAB, rhoCD, rhoACBD};
      size_t i = ndi_.to_linear(multi_index);
      return operator[](i);
    }

    double operator()(size_t rhoAC, size_t rhoBD, size_t rhoABCD, 
    				          size_t rhoAB, size_t rhoCD, size_t rhoACBD) const
    {
      std::array multi_index = {rhoAC, rhoBD, rhoABCD, rhoAB, rhoCD, rhoACBD};
      size_t i = ndi_.to_linear(multi_index);
      return operator[](i);
    }

    friend std::ostream& operator<<(std::ostream& s, const U3U9& x)
  	{
  	  s << "rhoAC\trhoBD\trhoABCD\trhoAB\trhoCD\trhoACBD\n";
  	  std::vector<size_t> multi(6);
  	  for (size_t i = 0; i < x.size(); i++)
  	  {
  	    x.ndi_.to_multi(i, multi);
  	    s << multi[0] << "\t" << multi[1] << "\t"
  	      << multi[2] << "\t" << multi[3] << "\t"
  	      << multi[4] << "\t" << multi[5] << "\t -> "
  	      << x[i] << "\n";
  	  }
	    return s;
  	}

    void print(std::ostream& s = std::cout) const
    {
      s << *this;
    }
  };

} // namespace u3

#endif
