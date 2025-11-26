#include "SpectralFft.h"
#include <cstring> // for memset

namespace samples {
namespace frostflow {

SpectralFft::SpectralFft() {
    setup = nullptr;
    workBuffer = nullptr;
}

SpectralFft::~SpectralFft() {
#if !defined(erb_TARGET_DAISY)
    if (setup) {
        pffft_destroy_setup(setup);
    }
    if (workBuffer) {
        pffft_aligned_free(workBuffer);
    }
#endif
}

void SpectralFft::init() {
#if !defined(erb_TARGET_DAISY)
    if (setup) return; // Already initialized
    setup = pffft_new_setup(FFT_SIZE, PFFFT_REAL);
    workBuffer = (float*)pffft_aligned_malloc(FFT_SIZE * sizeof(float));
#endif
}

void SpectralFft::forward(const float* input, float* output) {
#if !defined(erb_TARGET_DAISY)
    if (!setup) return;
    // Note: pffft_transform_ordered handles the ordering for us
    // input and output can point to the same buffer if needed, but we assume separate for clarity
    pffft_transform_ordered(setup, input, output, workBuffer, PFFFT_FORWARD);
#endif
}

void SpectralFft::inverse(const float* input, float* output) {
#if !defined(erb_TARGET_DAISY)
    if (!setup) return;
    pffft_transform_ordered(setup, input, output, workBuffer, PFFFT_BACKWARD);
    
    // PFFFT inverse transform is unscaled (result is N * original).
    // We need to scale by 1/N.
    const float scale = 1.0f / FFT_SIZE;
    for (size_t i = 0; i < FFT_SIZE; ++i) {
        output[i] *= scale;
    }
#endif
}

}
}
