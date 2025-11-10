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

inline std::vector<float> generateSawtooth(float frequencyHz, float amplitude, float sampleRate, size_t sampleCount)
{
    std::vector<float> buffer(sampleCount);
    const float period = sampleRate / frequencyHz;

    for (size_t i = 0; i < sampleCount; ++i) {
        const float t = static_cast<float>(i);
        const float phase = std::fmod(t, period) / period;
        buffer[i] = amplitude * (2.0f * phase - 1.0f);
    }

    return buffer;
}

inline std::vector<float> generateSineSweep(float startHz,
                                            float endHz,
                                            float amplitude,
                                            float sampleRate,
                                            size_t sampleCount)
{
    std::vector<float> buffer(sampleCount);
    float phase = 0.0f;

    for (size_t i = 0; i < sampleCount; ++i) {
        const float t = static_cast<float>(i) / static_cast<float>(sampleCount - 1);
        const float freq = startHz + (endHz - startHz) * t;
        const float phaseIncrement = 2.0f * kPi * freq / sampleRate;
        phase += phaseIncrement;
        buffer[i] = amplitude * std::sin(phase);
    }

    return buffer;
}

inline std::vector<float> generateMultiSine(const std::vector<float>& frequenciesHz,
                                            const std::vector<float>& amplitudes,
                                            float sampleRate,
                                            size_t sampleCount)
{
    std::vector<float> buffer(sampleCount, 0.0f);
    if (frequenciesHz.size() != amplitudes.size()) {
        return buffer;
    }

    std::vector<float> phases(frequenciesHz.size(), 0.0f);
    for (size_t i = 0; i < sampleCount; ++i) {
        float sample = 0.0f;
        for (size_t h = 0; h < frequenciesHz.size(); ++h) {
            const float phaseIncrement = 2.0f * kPi * frequenciesHz[h] / sampleRate;
            phases[h] += phaseIncrement;
            sample += amplitudes[h] * std::sin(phases[h]);
        }
        buffer[i] = sample;
    }

    return buffer;
}

inline std::vector<float> generateHarmonicSum(float fundamentalHz,
                                              const std::vector<std::size_t>& harmonics,
                                              float baseAmplitude,
                                              float sampleRate,
                                              size_t sampleCount)
{
    std::vector<float> buffer(sampleCount, 0.0f);
    if (harmonics.empty()) {
        return buffer;
    }

    for (std::size_t h : harmonics) {
        const float amplitude = baseAmplitude / static_cast<float>(harmonics.size());
        const float frequency = fundamentalHz * static_cast<float>(h);
        const auto partial = generateSine(frequency, amplitude, sampleRate, sampleCount);
        for (size_t i = 0; i < sampleCount; ++i) {
            buffer[i] += partial[i];
        }
    }

    return buffer;
}


