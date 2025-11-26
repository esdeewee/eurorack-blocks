#include "../test.h"
#include "TestHelpers.h"
#include "../../../samples/frostflow/dsp/SpectralFft.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <cstring>
#include <pffft.h>

using namespace samples::frostflow;

// Helper wrapper to use the TestHelpers allocators
struct FftBuffers {
    float* in;
    float* out;
    float* res;
    
    FftBuffers(size_t size) {
        in = (float*)pffft_aligned_malloc(size * sizeof(float));
        out = (float*)pffft_aligned_malloc(size * sizeof(float));
        res = (float*)pffft_aligned_malloc(size * sizeof(float));
        clear(size);
    }
    
    ~FftBuffers() {
        pffft_aligned_free(in);
        pffft_aligned_free(out);
        pffft_aligned_free(res);
    }
    
    void clear(size_t size) {
        memset(in, 0, size * sizeof(float));
        memset(out, 0, size * sizeof(float));
        memset(res, 0, size * sizeof(float));
    }
};

erb_TEST_CASE(FrostFlow, SpectralFft_Sine440) {
    SpectralFft fft;
    fft.init();
    FftBuffers buf(1024);
    
    TestSignal::generateSine(buf.in, 1024, 440.0f, 48000.0f);
    
    fft.forward(buf.in, buf.out);
    fft.inverse(buf.out, buf.res);
    
    float maxDiff = 0.0f;
    for (size_t i = 0; i < 1024; ++i) {
        float diff = std::abs(buf.in[i] - buf.res[i]);
        if (diff > maxDiff) maxDiff = diff;
    }
    
    erb_ASSERT_FLOAT_EQ(maxDiff, 0.0f, 0.001f, "Sine reconstruction failed");
}

erb_TEST_CASE(FrostFlow, SpectralFft_Impulse) {
    SpectralFft fft;
    fft.init();
    FftBuffers buf(1024);
    
    TestSignal::generateImpulse(buf.in, 1024, 1.0f);
    
    fft.forward(buf.in, buf.out);
    
    // Impulse in Time = Flat in Freq (Magnitude)
    // Check a few bins magnitudes (PFFFT ordered)
    // Bin 0: DC
    // Bin 1: Nyquist
    // Bin 2k, 2k+1: Real/Imag
    
    // For real FFT of [1, 0, 0...], all Real parts should be 1.0?
    // Let's check reconstruction first.
    
    fft.inverse(buf.out, buf.res);
    
    float maxDiff = 0.0f;
    for (size_t i = 0; i < 1024; ++i) {
        float diff = std::abs(buf.in[i] - buf.res[i]);
        if (diff > maxDiff) maxDiff = diff;
    }
    
    erb_ASSERT_FLOAT_EQ(maxDiff, 0.0f, 0.001f, "Impulse reconstruction failed");
}

erb_TEST_CASE(FrostFlow, SpectralFft_DC) {
    SpectralFft fft;
    fft.init();
    FftBuffers buf(1024);
    
    TestSignal::generateDC(buf.in, 1024, 0.5f);
    
    fft.forward(buf.in, buf.out);
    
    // DC Bin (Index 0) should be N * DC? Or Sum?
    // Sum of 0.5 * 1024 = 512.
    // PFFFT outputs sum.
    
    erb_ASSERT_FLOAT_EQ(buf.out[0], 512.0f, 0.1f, "DC Bin analysis failed");
    
    fft.inverse(buf.out, buf.res);
    
    float maxDiff = 0.0f;
    for (size_t i = 0; i < 1024; ++i) {
        float diff = std::abs(buf.in[i] - buf.res[i]);
        if (diff > maxDiff) maxDiff = diff;
    }
    
    erb_ASSERT_FLOAT_EQ(maxDiff, 0.0f, 0.001f, "DC reconstruction failed");
}

erb_TEST_CASE(FrostFlow, SpectralFft_Nyquist) {
    SpectralFft fft;
    fft.init();
    FftBuffers buf(1024);
    
    // Nyquist is [1, -1, 1, -1...] sequence
    for (size_t i=0; i<1024; ++i) buf.in[i] = (i % 2 == 0) ? 1.0f : -1.0f;
    
    fft.forward(buf.in, buf.out);
    
    // Nyquist Bin is Index 1 in PFFFT Ordered
    // Sum of absolute vals? 
    // It's basically cosine at Fs/2.
    // Should be high energy in bin 1.
    
    fft.inverse(buf.out, buf.res);
    
    float maxDiff = 0.0f;
    for (size_t i = 0; i < 1024; ++i) {
        float diff = std::abs(buf.in[i] - buf.res[i]);
        if (diff > maxDiff) maxDiff = diff;
    }
    
    erb_ASSERT_FLOAT_EQ(maxDiff, 0.0f, 0.001f, "Nyquist reconstruction failed");
}

erb_TEST_CASE(FrostFlow, SpectralFft_WhiteNoise) {
    SpectralFft fft;
    fft.init();
    FftBuffers buf(1024);
    
    TestSignal::generateNoise(buf.in, 1024);
    
    fft.forward(buf.in, buf.out);
    fft.inverse(buf.out, buf.res);
    
    float maxDiff = 0.0f;
    for (size_t i = 0; i < 1024; ++i) {
        float diff = std::abs(buf.in[i] - buf.res[i]);
        if (diff > maxDiff) maxDiff = diff;
    }
    
    erb_ASSERT_FLOAT_EQ(maxDiff, 0.0f, 0.001f, "Noise reconstruction failed");
}
