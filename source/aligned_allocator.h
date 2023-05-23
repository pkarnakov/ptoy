#pragma once

#include <cstdlib>

template <typename T, unsigned alignment>
class AlignedAllocator {
 public:
  using value_type = T;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  template <class U>
  struct rebind {
    using other = AlignedAllocator<U, alignment>;
  };
  pointer allocate(size_type n) {
    pointer p;
    if (posix_memalign(
            reinterpret_cast<void**>(&p), alignment, n * sizeof(T))) {
      throw std::bad_alloc();
    }
    return p;
  }
  void deallocate(pointer p, size_type) {
    std::free(p);
  }
  void construct(pointer p, const_reference val) {
    new (static_cast<void*>(p)) T(val);
  }
  void destroy(pointer p) {
    p->~T();
  }
  size_type max_size() const {
    std::allocator<T> a;
    return a.max_size();
  }
  pointer address(reference x) const;
  const_pointer address(const_reference x) const;

  bool operator==(const AlignedAllocator&) const {
    return true;
  }
  bool operator!=(const AlignedAllocator&) const {
    return false;
  }

  AlignedAllocator() {}
  AlignedAllocator(AlignedAllocator const&) {}

  template <typename U>
  AlignedAllocator(AlignedAllocator<U, alignment> const&) noexcept {}
};
