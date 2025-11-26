#pragma once
#include <vector>
#include <memory>
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>
#include <pffft.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

// --- Memory Helpers ---

// Simple Aligned Allocator for std::vector
template <typename T, std::size_t Alignment = 32>
struct AlignedAllocator {
    using value_type = T;

    AlignedAllocator() = default;
    template <class U> AlignedAllocator(const AlignedAllocator<U, Alignment>&) {}

    T* allocate(std::size_t n) {
        return (T*)pffft_aligned_malloc(n * sizeof(T));
    }

    void deallocate(T* p, std::size_t) {
        pffft_aligned_free(p);
    }
};

using AlignedFloatVector = std::vector<float, AlignedAllocator<float, 32>>;


// --- Signal Generators ---

struct TestSignal {
    static void generateSine(float* buffer, size_t size, float freq, float sampleRate, float amplitude = 1.0f) {
        for (size_t i = 0; i < size; ++i) {
            buffer[i] = amplitude * std::sin(2.0f * M_PI * freq * i / sampleRate);
        }
    }

    static void generateCosine(float* buffer, size_t size, float freq, float sampleRate, float amplitude = 1.0f) {
        for (size_t i = 0; i < size; ++i) {
            buffer[i] = amplitude * std::cos(2.0f * M_PI * freq * i / sampleRate);
        }
    }
    
    static void generateDC(float* buffer, size_t size, float level) {
        for (size_t i = 0; i < size; ++i) {
            buffer[i] = level;
        }
    }
    
    static void generateImpulse(float* buffer, size_t size, float amplitude = 1.0f) {
        std::fill(buffer, buffer + size, 0.0f);
        if (size > 0) buffer[0] = amplitude;
    }

    static void generateNoise(float* buffer, size_t size, float amplitude = 1.0f, int seed = 42) {
        srand(seed);
        for (size_t i = 0; i < size; ++i) {
            buffer[i] = amplitude * ((float)rand() / RAND_MAX * 2.0f - 1.0f);
        }
    }

    static void generateSaw(float* buffer, size_t size, float freq, float sampleRate, float amplitude = 1.0f) {
         for (size_t i = 0; i < size; ++i) {
            float phase = std::fmod((float)i * freq / sampleRate, 1.0f);
            buffer[i] = amplitude * (2.0f * phase - 1.0f);
        }
    }
};

// --- Assertions with Logging ---

#define erb_ASSERT_MSG(expr, msg) \
    do { \
        if (!(expr)) { \
            std::cout << "[ASSERT FAILURE] " << msg << std::endl; \
            std::cout << "      File: " << __FILE__ << ", Line: " << __LINE__ << std::endl; \
            std::terminate(); \
        } \
    } while(false)

#define erb_ASSERT_FLOAT_EQ(val1, val2, tolerance, context) \
    do { \
        float diff = std::abs((val1) - (val2)); \
        if (diff > (tolerance)) { \
            std::cout << "[FLOAT MISMATCH] " << context << std::endl; \
            std::cout << "      Expected: " << (val2) << ", Got: " << (val1) << ", Diff: " << diff << std::endl; \
            std::terminate(); \
        } \
    } while(false)
