#ifndef UTILS_OPTIONAL_REFERENCE_H
#define UTILS_OPTIONAL_REFERENCE_H

#include <cassert>

namespace utils
{
  // reassignable reference wrapper with empty state:
  template <typename T>
  class optional_reference
  {
    T* ptr_;

  public:
    optional_reference() : ptr_(nullptr) { }

    optional_reference(const optional_reference&) = default;

    template <typename U>
    optional_reference(U& obj) : ptr_(&obj)
    { }

    optional_reference& operator=(const optional_reference&) = default;

    template <typename U>
    optional_reference& operator=(U& obj)
    {
      ptr_ = &obj;
      return *this;
    }

    bool empty() const { return ptr_ == nullptr; }

    T& get() const
    {
      assert(!empty());
      return *ptr_;
    }

    operator T&() const { return get(); }

    T* operator->() const
    {
      assert(!empty());
      return ptr_;
    }

    void reset() { ptr_ = nullptr; }
  };

  template <typename T>
  inline auto ref(T& r)
  {
    return optional_reference<T>(r);
  }

  template <typename T>
  inline auto cref(const T& r)
  {
    return optional_reference<const T>(r);
  }

} // namespace utils

#endif
