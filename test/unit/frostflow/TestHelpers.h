#pragma once
#include <vector>
#include <memory>
#include <cstdlib>

// Simple Aligned Allocator for std::vector
template <typename T, std::size_t Alignment = 16>
struct AlignedAllocator {
    using value_type = T;

    AlignedAllocator() = default;
    template <class U> AlignedAllocator(const AlignedAllocator<U, Alignment>&) {}

    T* allocate(std::size_t n) {
        void* ptr = nullptr;
#if defined(_MSC_VER)
        ptr = _aligned_malloc(n * sizeof(T), Alignment);
#else
        if (posix_memalign(&ptr, Alignment, n * sizeof(T)) != 0) ptr = nullptr;
#endif
        if (!ptr) throw std::bad_alloc();
        return static_cast<T*>(ptr);
    }

    void deallocate(T* p, std::size_t) {
#if defined(_MSC_VER)
        _aligned_free(p);
#else
        free(p);
#endif
    }
};

using AlignedFloatVector = std::vector<float, AlignedAllocator<float, 16>>;

