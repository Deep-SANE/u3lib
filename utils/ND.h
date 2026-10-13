#ifndef UTILS_ND_H
#define UTILS_ND_H

#include <cassert>
#include <cstddef>
#include <limits>
#include <span>
#include <type_traits>
#include <vector>

#include <boost/container/static_vector.hpp>

namespace utils
{
  // template argument for dynamic number of dimensions:
  inline constexpr size_t DYNAMIC = std::numeric_limits<size_t>::max();

  template <size_t NDIM = DYNAMIC>
    requires(NDIM > 0)
  class NDIndexT
  {
    // type for storing dimensions:
    using vector_type =
        std::conditional_t<NDIM == DYNAMIC, std::vector<size_t>,
                           boost::container::static_vector<size_t, NDIM>>;

    vector_type dimensions_;
    size_t cardinality_ = 0;

  public:
    size_t n_dimensions() const { return dimensions_.size(); }

    size_t dimension(size_t i) const
    {
      assert(i < n_dimensions());
      return dimensions_[i];
    }

    size_t cardinality() const { return cardinality_; }

    void add_dimension(size_t N)
    {
      assert(n_dimensions() < NDIM);

      if (n_dimensions() == 0)
        cardinality_ = N;
      else
        cardinality_ *= N;

      dimensions_.push_back(N);
    }

    void add_dimensions(auto N) { add_dimension(N); }

    void add_dimensions(auto head, auto... tail)
    {
      add_dimension(head);
      add_dimensions(tail...);
    }

    template <typename T>
    size_t to_linear(const T& multi_index) const
    {
      assert(cardinality() > 0);
      assert(n_dimensions() > 0);
      assert(n_dimensions() == multi_index.size());

      size_t linear_index = 0;
      size_t stride = 1;

      for (long i = n_dimensions() - 1; i >= 0; i--)
      {
        assert(dimension(i) > 0);

        linear_index += multi_index[i] * stride;
        stride *= dimension(i);
      }

      return linear_index;
    }

    // result is passed as a span by value, which allows to reuse memory
    template <typename T>
    void to_multi(size_t linear_index, T& multi_index) const
    {
      assert(cardinality() > 0);
      assert(n_dimensions() > 0);
      assert(n_dimensions() == multi_index.size());

      for (long i = n_dimensions() - 1; i >= 0; i--)
      {
        assert(dimension(i) > 0);

        multi_index[i] = linear_index % dimension(i);
        linear_index /= dimension(i);
      }
    }
  };

  using NDIndex = NDIndexT<>;
  using NDIndex2 = NDIndexT<2>;
  using NDIndex3 = NDIndexT<3>;
  using NDIndex4 = NDIndexT<4>;
  using NDIndex5 = NDIndexT<5>;
  using NDIndex6 = NDIndexT<6>;
  using NDIndex7 = NDIndexT<7>;
  using NDIndex8 = NDIndexT<8>;
  using NDIndex9 = NDIndexT<9>;
  using NDIndex10 = NDIndexT<10>;
  using NDIndex11 = NDIndexT<11>;
  using NDIndex12 = NDIndexT<12>;
  using NDIndex13 = NDIndexT<13>;
  using NDIndex14 = NDIndexT<14>;
  using NDIndex15 = NDIndexT<15>;
  using NDIndex16 = NDIndexT<16>;
  using NDIndex17 = NDIndexT<17>;
  using NDIndex18 = NDIndexT<18>;
  using NDIndex19 = NDIndexT<19>;
  using NDIndex20 = NDIndexT<20>;

  template <typename T>
  class NDSpan
  {
    NDIndex ndi_;
    std::vector<std::span<const T>> spans_;

    std::vector<size_t> multi_index_;

  public:
    template <typename U>
    void add_dimension(const U& elements)
    {
      spans_.emplace_back(elements);
      multi_index_.emplace_back();
      ndi_.add_dimension(spans_.back().size());
    }

    size_t n_dimensions() const { return spans_.size(); }

    size_t cardinality() const { return ndi_.cardinality(); }

    template <typename U>
    void extract_multi_element(size_t linear_index, U& where)
    {
      assert(std::size(where) == n_dimensions());
      assert(multi_index_.size() == n_dimensions());

      ndi_.to_multi(linear_index, multi_index_);

      for (int i = 0; i < n_dimensions(); i++)
      {
        auto index = multi_index_[i];
        assert(index < spans_[i].size());

        where[i] = spans_[i][index];
      }
    }
  };

  template <typename T, std::size_t NDIM>
  struct NDVector
  {
    std::vector<T> vector;
    NDIndexT<NDIM> index;
  };

} // namespace utils

#endif
