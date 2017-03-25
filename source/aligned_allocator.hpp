// Example codes for HPC course
// (c) 2012/15 Andreas Hehn, ETH Zurich

#ifndef HPC12_ALIGNED_ALLOCATOR_HPP
#define HPC12_ALIGNED_ALLOCATOR_HPP

#ifdef _WIN32
#include <malloc.h>
#else
#include <cstdlib>
#endif

#if __cplusplus >= 201103L
#define NOEXCEPT_SPEC noexcept
#else
#define NOEXCEPT_SPEC
#endif

namespace hpc15 {

// Alignment must be a power of 2 !
template <typename T, unsigned int Alignment>
class aligned_allocator {
  public:
    typedef T*              pointer;
    typedef T const*        const_pointer;
    typedef T&              reference;
    typedef T const&        const_reference;
    typedef T               value_type;
    typedef std::size_t     size_type;
    typedef std::ptrdiff_t  difference_type;

    template <typename U>
    struct rebind {
        typedef aligned_allocator<U,Alignment> other;
    };

    aligned_allocator() NOEXCEPT_SPEC {
    }

    aligned_allocator(aligned_allocator const& a) NOEXCEPT_SPEC {
    }

    template <typename U>
    aligned_allocator(aligned_allocator<U,Alignment> const& b) NOEXCEPT_SPEC {
    }

    pointer allocate(size_type n) {
        pointer p;
#ifdef _WIN32
        p = _aligned_malloc(n*sizeof(T), Alignment);
        if(p == 0)
            throw std::bad_alloc();
#else
        if(posix_memalign(reinterpret_cast<void**>(&p), Alignment, n * sizeof(T) ))
            throw std::bad_alloc();
#endif
        return p;
    }

    void deallocate(pointer p, size_type n) NOEXCEPT_SPEC {
        std::free(p);
    }

    size_type max_size() const NOEXCEPT_SPEC {
        std::allocator<T> a;
        return a.max_size();
    }

#if __cplusplus >= 201103L
    template <typename C, class... Args>
    void construct(C* c, Args&&... args) {
        new ((void*)c) C(std::forward<Args>(args)...);
    }
#else
    void construct(pointer p, const_reference t) {
        new((void *)p) T(t);
    }
#endif

    template <typename C>
    void destroy(C* c) {
        c->~C();
    }

    bool operator == (aligned_allocator const& a2) const NOEXCEPT_SPEC {
        return true;
    }

    bool operator != (aligned_allocator const& a2) const NOEXCEPT_SPEC {
        return false;
    }

    template <typename U, unsigned int UAlignment>
    bool operator == (aligned_allocator<U,UAlignment> const& b) const NOEXCEPT_SPEC {
        return false;
    }

    template <typename U, unsigned int UAlignment>
    bool operator != (aligned_allocator<U,UAlignment> const& b) const NOEXCEPT_SPEC {
        return true;
    }
};

}

#undef NOEXPECT_SPEC

#endif //HPC12_ALIGNED_ALLOCATOR_HPP
