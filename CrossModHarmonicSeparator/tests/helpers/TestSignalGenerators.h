#pragma once

#include <vector>
#include <cmath>

constexpr float kPi = 3.14159265358979323846f;

inline std::vector<float> generateSine(float frequencyHz, float amplitude, float sampleRate, size_t sampleCount)
{
    std::vector<float> buffer(sampleCount);
    const float phaseIncrement = 2.0f * kPi * frequencyHz / sampleRate;

    float phase = 0.0f;
    for (size_t i = 0; i < sampleCount; ++i) {
        buffer[i] = amplitude * std::sin(phase);
        phase += phaseIncrement;
    }

    return buffer;
}


