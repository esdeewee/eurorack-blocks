#include "SpectralFreeze.h"
#include <cmath>
#include <cstring>

#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif

namespace samples {
namespace frostflow {

SpectralFreeze::SpectralFreeze() 
    : hopCounter(0), wasFrozen(false), fftWorkIn(nullptr), fftWorkOut(nullptr)
{
    inputBufferL.resize(FFT_SIZE, 0.0f);
    inputBufferR.resize(FFT_SIZE, 0.0f);
    outputBufferL.resize(FFT_SIZE, 0.0f);
    outputBufferR.resize(FFT_SIZE, 0.0f);
    
    currentMagL.resize(FFT_SIZE, 0.0f);
    currentPhaseL.resize(FFT_SIZE, 0.0f);
    currentMagR.resize(FFT_SIZE, 0.0f);
    currentPhaseR.resize(FFT_SIZE, 0.0f);
    
    frozenMagL.resize(FFT_SIZE, 0.0f);
    frozenPhaseL.resize(FFT_SIZE, 0.0f);
    frozenMagR.resize(FFT_SIZE, 0.0f);
    frozenPhaseR.resize(FFT_SIZE, 0.0f);

#if !defined(erb_TARGET_DAISY)
    fftWorkIn = (float*)pffft_aligned_malloc(FFT_SIZE * sizeof(float));
    fftWorkOut = (float*)pffft_aligned_malloc(FFT_SIZE * sizeof(float));
#else
    // Daisy: Use SRAM or simple malloc if handled elsewhere
    // Phase 2: We will handle this.
#endif
    
    // Initialize Hann Window
    window.resize(FFT_SIZE);
    for (size_t i = 0; i < FFT_SIZE; ++i) {
        window[i] = 0.5f * (1.0f - cosf(2.0f * M_PI * i / (float)(FFT_SIZE - 1)));
    }
}

SpectralFreeze::~SpectralFreeze() {
#if !defined(erb_TARGET_DAISY)
    if (fftWorkIn) pffft_aligned_free(fftWorkIn);
    if (fftWorkOut) pffft_aligned_free(fftWorkOut);
#endif
}

void SpectralFreeze::init() {
    fft.init();
    processor.init();
}

void SpectralFreeze::process(float inL, float inR, float& outL, float& outR, bool freeze) {
    // DUMMY COMMENT 1
    // DUMMY COMMENT 2
    // DUMMY COMMENT 3
    // DUMMY COMMENT 4
    // DUMMY COMMENT 5
    // DUMMY COMMENT 6
    // DUMMY COMMENT 7
    // DUMMY COMMENT 8
    // DUMMY COMMENT 9
    // DUMMY COMMENT 10
    // 1. Push inputs to FIFO (Shift and Push)
    
    // Read output current sample
    outL = outputBufferL[0];
    outR = outputBufferR[0];

    // Shift output buffer
    memmove(outputBufferL.data(), outputBufferL.data() + 1, (FFT_SIZE - 1) * sizeof(float));
    memmove(outputBufferR.data(), outputBufferR.data() + 1, (FFT_SIZE - 1) * sizeof(float));
    outputBufferL[FFT_SIZE - 1] = 0.0f;
    outputBufferR[FFT_SIZE - 1] = 0.0f;

    // Shift input buffer and push new sample
    memmove(inputBufferL.data(), inputBufferL.data() + 1, (FFT_SIZE - 1) * sizeof(float));
    memmove(inputBufferR.data(), inputBufferR.data() + 1, (FFT_SIZE - 1) * sizeof(float));
    inputBufferL[FFT_SIZE - 1] = inL;
    inputBufferR[FFT_SIZE - 1] = inR;

    hopCounter++;
    if (hopCounter >= HOP_SIZE) {
        hopCounter = 0;
        
        // Handle Freeze State Logic
        if (freeze && !wasFrozen) {
            // Just entered freeze state: Capture current spectrum
            // We need to analyze the current buffer once to get the freeze frame
            analyze(inputBufferL, frozenMagL, frozenPhaseL);
            analyze(inputBufferR, frozenMagR, frozenPhaseR);
        }
        wasFrozen = freeze;

        if (freeze) {
            // Synthesis from Frozen Data
            // Apply processor (motion/decay) to a copy of the frozen data
            // Copy frozen data to current work buffers so processor can modify them non-destructively
            currentMagL = frozenMagL;
            currentPhaseL = frozenPhaseL;
            currentMagR = frozenMagR;
            currentPhaseR = frozenPhaseR;
            
            // Apply Effects
            processor.process(currentMagL.data(), currentPhaseL.data(), FFT_SIZE);
            processor.process(currentMagR.data(), currentPhaseR.data(), FFT_SIZE);
            
            // Resynthesize
            synthesize(outputBufferL, currentMagL, currentPhaseL);
            synthesize(outputBufferR, currentMagR, currentPhaseR);
        } else {
            // Live Processing
            // Analyze -> Resynthesize (Passthrough with spectral effects)
            analyze(inputBufferL, currentMagL, currentPhaseL);
            analyze(inputBufferR, currentMagR, currentPhaseR);
            
            // Apply Effects even in live mode? Or bypass?
            // Specification says: "Works on frozen layers only (doesn't affect live input)"
            // So in live mode, we bypass processor.
            
            synthesize(outputBufferL, currentMagL, currentPhaseL);
            synthesize(outputBufferR, currentMagR, currentPhaseR);
        }
    }
}

void SpectralFreeze::analyze(const std::vector<float>& input, std::vector<float>& mag, std::vector<float>& phase) {
    // Apply Window
    for (size_t i = 0; i < FFT_SIZE; ++i) {
        fftWorkIn[i] = input[i] * window[i];
    }

    // FFT
    fft.forward(fftWorkIn, fftWorkOut);

    // Cartesian to Polar
    mag[0] = std::abs(fftWorkOut[0]);
    phase[0] = 0.0f; // DC phase 0

    mag[1] = std::abs(fftWorkOut[1]);
    phase[1] = 0.0f; // Nyquist phase 0

    for (size_t k = 1; k < FFT_SIZE / 2; ++k) {
        float re = fftWorkOut[2 * k];
        float im = fftWorkOut[2 * k + 1];
        mag[2 * k] = std::sqrt(re * re + im * im);
        phase[2 * k] = std::atan2(im, re);
        
        mag[2 * k + 1] = mag[2 * k];
        phase[2 * k + 1] = phase[2 * k];
    }
}

void SpectralFreeze::synthesize(std::vector<float>& output, const std::vector<float>& mag, const std::vector<float>& phase) {
    // Polar to Cartesian
    fftWorkIn[0] = mag[0]; // DC Real only
    fftWorkIn[1] = mag[1]; // Nyquist Real only

    for (size_t k = 1; k < FFT_SIZE / 2; ++k) {
        float m = mag[2 * k];
        float p = phase[2 * k];
        
        fftWorkIn[2 * k] = m * std::cos(p);
        fftWorkIn[2 * k + 1] = m * std::sin(p);
    }

    // Inverse FFT
    fft.inverse(fftWorkIn, fftWorkOut);

    // Window and Overlap-Add
    for (size_t i = 0; i < FFT_SIZE; ++i) {
        output[i] += fftWorkOut[i] * window[i] * (2.0f / 3.0f);
    }
}

}
}
