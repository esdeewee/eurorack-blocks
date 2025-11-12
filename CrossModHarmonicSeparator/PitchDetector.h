#pragma once

#include <cstddef>
#include <optional>
#include <vector>

class PitchDetector
{
public:
    PitchDetector(float sampleRate, std::size_t analysisSize = 2048, std::size_t hopSize = 48);

    void reset();

    void processBuffer(const float* buffer, std::size_t length);

    bool hasPitch() const;
    float getCurrentPitchHz() const;
    float getLastPitchHz() const;

private:
    bool detectPitch(const float* windowStart, std::size_t windowLength, float& outPitchHz);
    float computeRms(const float* windowStart, std::size_t length) const;

    float sample_rate_;
    std::size_t analysis_size_;
    std::size_t hop_size_;

    std::vector<float> accumulator_;
    std::vector<float> difference_;
    std::vector<float> cmnd_;

    std::optional<float> current_pitch_hz_;
    std::optional<float> last_valid_pitch_hz_;
};


