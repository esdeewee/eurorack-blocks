#pragma once

#include <cstdint>

namespace samples {
namespace frostflow {

class RngWrapper {
public:
    RngWrapper() : state(123456789) {}

    void setSeed(uint32_t seed) {
        if (seed == 0) seed = 123456789;
        state = seed;
    }

    // Returns float in [0, 1)
    float nextFloat() {
        uint32_t x = nextInt();
        // Convert to float: (x >> 9) | 0x3f800000 gives float in [1, 2)
        // Then subtract 1.0f.
        // This is fast standard trick.
        union {
            uint32_t i;
            float f;
        } u;
        u.i = (x >> 9) | 0x3f800000;
        return u.f - 1.0f;
    }
    
    // Returns float in [-1, 1)
    float nextBiFloat() {
        return nextFloat() * 2.0f - 1.0f;
    }

private:
    uint32_t nextInt() {
        // Xorshift32
        uint32_t x = state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        state = x;
        return x;
    }

    uint32_t state;
};

}
}
