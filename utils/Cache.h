#ifndef UTILS_CACHE_H
#define UTILS_CACHE_H

#include <cassert>
#include <utility>

#ifdef UTILS_CACHE_LOG
#include <map>
#include <source_location>
#endif

#include <boost/container_hash/hash.hpp>

#ifdef USE_OPENMP
#define UTILS_CACHE_USE_CONCURRENT_FLATMAP
#endif

//#define UTILS_CACHE_USE_FLATMAP

#if defined  UTILS_CACHE_USE_CONCURRENT_FLATMAP
#include <boost/unordered/concurrent_flat_map.hpp>
#elif defined UTILS_CACHE_USE_FLATMAP
#include <boost/unordered/unordered_flat_map.hpp>
#else
#include <unordered_map>
#endif

namespace utils
{
  // Class that implements caching/memoization based on a hash table:
  template <typename Key, typename Value, typename Hash = boost::hash<Key>>
  class Cache
  {
#if defined UTILS_CACHE_USE_CONCURRENT_FLATMAP
    using map_type = boost::concurrent_flat_map<Key, Value, Hash>;
#elif defined UTILS_CACHE_USE_FLATMAP
    using map_type = boost::unordered_flat_map<Key, Value, Hash>;
#else
    using map_type = std::unordered_map<Key, Value, Hash>;
#endif

    map_type map_;

#ifdef UTILS_CACHE_LOG
    struct helper
    {
      std::map<std::string, std::size_t> m;

      ~helper()
      {
        for (const auto& [func, count] : m)
          std::cerr << func << " : " << count << "\n";
      }
    };
#endif

  public:
    void clear()
    {
      map_type().swap(map_);
    }
    
    template <typename Prod, typename Acc>
    void visit(const Key& key, Prod&& prod, Acc&& acc
#ifdef UTILS_CACHE_LOG
               ,
               std::source_location sl = std::source_location::current()
#endif
    )
    {
#ifdef UTILS_CACHE_LOG
      static helper h;
      h.m[sl.function_name()]++;
#endif

#if defined UTILS_CACHE_USE_CONCURRENT_FLATMAP
      auto visited = map_.visit(key, std::forward<Acc>(acc));

      if (visited == 0)
      {
        auto value = std::forward<Prod>(prod)(key);

        map_.emplace(key, std::move(value));
        map_.visit(key, std::forward<Acc>(acc));
      }
#elif defined UTILS_CACHE_USE_FLATMAP
      auto iter = map_.find(key);
      if (iter == map_.end())
      {
        auto [iter2, inserted] =
            map_.try_emplace(key, std::forward<Prod>(prod)(key));

        assert(inserted);

        iter = iter2;
      }

      std::forward<Acc>(acc)(*iter);
#else
      auto iter = map_.find(key);
      if (iter == map_.end())
      {
        auto [iter2, inserted] =
            map_.try_emplace(key, std::forward<Prod>(prod)(key));

        assert(inserted);

        iter = iter2;
      }

      std::forward<Acc>(acc)(*iter);
#endif
    }
  };

} // namespace utils

#endif
