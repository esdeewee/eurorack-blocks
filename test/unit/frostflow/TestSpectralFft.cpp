#include "../test.h"
#include "../../../samples/frostflow/dsp/SpectralFft.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <cstring>
#include <pffft.h>
#include <cstdlib>

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, SpectralFft_RoundTrip) {
    SpectralFft fft;
    fft.init();
    
    float* input = (float*)pffft_aligned_malloc(1024 * sizeof(float));
    float* output = (float*)pffft_aligned_malloc(1024 * sizeof(float));
    float* result = (float*)pffft_aligned_malloc(1024 * sizeof(float));
    
    memset(input, 0, 1024 * sizeof(float));
    memset(output, 0, 1024 * sizeof(float));
    memset(result, 0, 1024 * sizeof(float));
    
    for (size_t i = 0; i < 1024; ++i) {
        input[i] = std::sin(2.0f * 3.14159f * 440.0f * i / 48000.0f);
    }
    
    fft.forward(input, output);
    
    if (std::abs(output[0]) > 1e10 || std::isnan(output[0])) std::cout << "Output[0] garbage: " << output[0] << std::endl;
    
    fft.inverse(output, result);
    
    if (std::abs(result[0]) > 1e10 || std::isnan(result[0])) std::cout << "Result[0] garbage: " << result[0] << std::endl;
    
    float maxDiff = 0.0f;
    for (size_t i = 0; i < 1024; ++i) {
        float diff = std::abs(input[i] - result[i]);
        if (diff > maxDiff) maxDiff = diff;
    }
    
    pffft_aligned_free(input);
    pffft_aligned_free(output);
    pffft_aligned_free(result);
    
    if (maxDiff >= 0.001f) {
        std::cout << "Max Diff: " << maxDiff << std::endl;
    }
    
    // Floating point error accumulation
    erb_TEST(maxDiff < 0.001f);
}

erb_TEST_CASE(FrostFlow, SpectralFft_WhiteNoise) {
    SpectralFft fft;
    fft.init();
    
    float* input = (float*)pffft_aligned_malloc(1024 * sizeof(float));
    float* output = (float*)pffft_aligned_malloc(1024 * sizeof(float));
    float* result = (float*)pffft_aligned_malloc(1024 * sizeof(float));
    
    // Deterministic seed
    srand(123);
    
    for (size_t i = 0; i < 1024; ++i) {
        input[i] = (float)rand() / RAND_MAX * 2.0f - 1.0f;
    }
    
    fft.forward(input, output);
    fft.inverse(output, result);
    
    float maxDiff = 0.0f;
    for (size_t i = 0; i < 1024; ++i) {
        float diff = std::abs(input[i] - result[i]);
        if (diff > maxDiff) maxDiff = diff;
    }
    
    pffft_aligned_free(input);
    pffft_aligned_free(output);
    pffft_aligned_free(result);
    
    if (maxDiff >= 0.001f) {
        std::cout << "Noise Max Diff: " << maxDiff << std::endl;
    }
    
    erb_TEST(maxDiff < 0.001f);
}
