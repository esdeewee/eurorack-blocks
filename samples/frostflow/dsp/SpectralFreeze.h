#pragma once

#include "SpectralFft.h"
#include "SpectralProcessor.h"
#include <vector>
#include <cmath>

namespace samples {
namespace frostflow {

class SpectralFreeze {
public:
    static const size_t FFT_SIZE = 1024;
    static const size_t HOP_SIZE = 256; // 4x overlap

    SpectralFreeze();
    ~SpectralFreeze();

    void init();

    /**
     * Process a block of audio.
     * @param inL Single sample left input
     * @param inR Single sample right input
     * @param outL Single sample left output
     * @param outR Single sample right output
     * @param freeze Gate state for freeze (true = frozen)
     */
    void process(float inL, float inR, float& outL, float& outR, bool freeze);

    // Parameter setters
    void setTilt(float tilt) { processor.setTilt(tilt); }
    void setDecay(float decay) { processor.setDecay(decay); }
    
private:
    void processFrame();
    void analyze(const std::vector<float>& input, std::vector<float>& mag, std::vector<float>& phase);
    void synthesize(std::vector<float>& output, const std::vector<float>& mag, const std::vector<float>& phase);

    SpectralFft fft;
    SpectralProcessor processor;

    // Buffers
    std::vector<float> inputBufferL;
    std::vector<float> inputBufferR;
    std::vector<float> outputBufferL;
    std::vector<float> outputBufferR;
    std::vector<float> window;

    // Spectral Data (Magnitudes and Phases)
    // We store separate Left and Right
    std::vector<float> currentMagL;
    std::vector<float> currentPhaseL;
    std::vector<float> currentMagR;
    std::vector<float> currentPhaseR;

    // Frozen Data
    std::vector<float> frozenMagL;
    std::vector<float> frozenPhaseL; // Initial phase at freeze point
    std::vector<float> frozenMagR;
    std::vector<float> frozenPhaseR;

    // Work buffers for FFT (Must be aligned)
    float* fftWorkIn;
    float* fftWorkOut;

    size_t hopCounter;
    bool wasFrozen;
};

}
}
