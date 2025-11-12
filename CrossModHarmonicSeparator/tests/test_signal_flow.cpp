#include <gtest/gtest.h>

#define CROSSMOD_TEST_ENV_HEADER "CrossModTestEnv.h"

#include <array>
#include <vector>
#include <numeric>

#include "CrossModHarmonicSeparator.h"
#include "helpers/DSPTestUtils.h"
#include "helpers/TestSignalGenerators.h"

namespace
{

constexpr float kSampleRate = 48000.0f;

} // namespace

TEST(SignalFlow, PassthroughSine)
{
    CrossModHarmonicSeparator module;
    module.init ();

    auto input = generateSine (1000.0f, 0.25f, kSampleRate, erb_BUFFER_SIZE);

    for (std::size_t i = 0; i < input.size (); ++i)
    {
        module.ui.audio_in [i] = input [i];
    }

    module.process ();

    std::vector<float> output (erb_BUFFER_SIZE);
    for (std::size_t i = 0; i < output.size (); ++i)
    {
        output [i] = module.ui.audio_out [i];
    }

    // Correlation should be near 1.
    double dot = 0.0;
    double inEnergy = 0.0;
    double outEnergy = 0.0;
    for (std::size_t i = 0; i < input.size (); ++i)
    {
        dot += static_cast<double>(input [i]) * static_cast<double>(output [i]);
        inEnergy += static_cast<double>(input [i]) * static_cast<double>(input [i]);
        outEnergy += static_cast<double>(output [i]) * static_cast<double>(output [i]);
    }

    double correlation = dot / std::sqrt (inEnergy * outEnergy);
    EXPECT_NEAR (correlation, 1.0, 1e-6);

    // RMS difference should be extremely small.
    std::vector<float> diff (erb_BUFFER_SIZE);
    for (std::size_t i = 0; i < diff.size (); ++i)
    {
        diff [i] = output [i] - input [i];
    }
    EXPECT_LT (computeRms (diff), 1e-6);
}

TEST(SignalFlow, ImpulseLatency)
{
    CrossModHarmonicSeparator module;
    module.init ();

    std::vector<float> impulse (erb_BUFFER_SIZE, 0.0f);
    impulse [0] = 1.0f;

    for (std::size_t i = 0; i < impulse.size (); ++i)
    {
        module.ui.audio_in [i] = impulse [i];
    }

    module.process ();

    std::size_t firstOutputIndex = erb_BUFFER_SIZE; // default to "not found"
    for (std::size_t i = 0; i < impulse.size (); ++i)
    {
        if (std::abs (module.ui.audio_out [i]) > 0.5f)
        {
            firstOutputIndex = i;
            break;
        }
    }

    // For direct passthrough we expect zero additional latency.
    EXPECT_EQ (firstOutputIndex, 0u);
}

TEST(SignalFlow, LowFrequencySeparationEnergy)
{
    CrossModHarmonicSeparator module;
    module.init();

    const float frequency = 60.0f;
    auto source = generateSine(frequency, 0.5f, kSampleRate, 48000);

    std::vector<float> routedInput;
    std::vector<float> separatedSum;

    std::size_t cursor = 0;
    float lastDetectedPitch = 0.0f;
    while (cursor < source.size()) {
        for (std::size_t i = 0; i < erb_BUFFER_SIZE; ++i) {
            const float sample = (cursor < source.size()) ? source[cursor++] : 0.0f;
            module.ui.audio_in[i] = sample;
            routedInput.push_back(sample);
        }

        module.process();
        if (module.has_detected_pitch) {
            lastDetectedPitch = module.detected_pitch_hz;
        }

        for (std::size_t i = 0; i < erb_BUFFER_SIZE; ++i) {
            const float separated = module.ui.fundamental_debug[i]
                                    + module.ui.prime_debug[i]
                                    + module.ui.composite_debug[i];
            separatedSum.push_back(separated);
        }
    }

    EXPECT_TRUE(module.has_spectral_separation);

    const std::size_t analysis = module.spectral_separator.analysisSize();
    const std::size_t hop = module.spectral_separator.hopSize();
    const std::size_t delay = (analysis > hop) ? (analysis - hop) : 0;

    ASSERT_GT(routedInput.size(), delay + 1024);
    ASSERT_GT(separatedSum.size(), delay + 1024);

    const std::size_t length = std::min(routedInput.size() - delay, separatedSum.size());
    double energyFundamental = 0.0;
    double energySum = 0.0;
    for (std::size_t i = delay; i < delay + length; ++i) {
        energyFundamental += static_cast<double>(routedInput[i]) * static_cast<double>(routedInput[i]);
    }
    for (std::size_t i = 0; i < length; ++i) {
        energySum += static_cast<double>(separatedSum[i]) * static_cast<double>(separatedSum[i]);
    }

    EXPECT_NEAR(energySum, energyFundamental, energyFundamental * 0.02);
}

TEST(SignalFlow, HighFrequencyReconstructionMatchesInput)
{
    CrossModHarmonicSeparator module;
    module.init();

    const float frequency = 1530.0f;
    auto source = generateSine(frequency, 0.4f, kSampleRate, 48000);

    std::vector<float> routedInput;
    std::vector<float> separatedSum;

    std::size_t cursor = 0;
    float lastDetectedPitch = 0.0f;
    while (cursor < source.size()) {
        for (std::size_t i = 0; i < erb_BUFFER_SIZE; ++i) {
            const float sample = (cursor < source.size()) ? source[cursor++] : 0.0f;
            module.ui.audio_in[i] = sample;
            routedInput.push_back(sample);
        }

        module.process();
        if (module.has_detected_pitch) {
            lastDetectedPitch = module.detected_pitch_hz;
        }

        for (std::size_t i = 0; i < erb_BUFFER_SIZE; ++i) {
            const float separated = module.ui.fundamental_debug[i]
                                    + module.ui.prime_debug[i]
                                    + module.ui.composite_debug[i];
            separatedSum.push_back(separated);
        }
    }

    EXPECT_TRUE(module.has_spectral_separation);

    const std::size_t analysis = module.spectral_separator.analysisSize();
    const std::size_t hop = module.spectral_separator.hopSize();
    const std::size_t delay = (analysis > hop) ? (analysis - hop) : 0;

    ASSERT_GT(routedInput.size(), delay + 1024);
    ASSERT_GT(separatedSum.size(), delay + 1024);

    const std::size_t length = std::min(routedInput.size() - delay, separatedSum.size());
    double diffEnergy = 0.0;
    double referenceEnergy = 0.0;
    for (std::size_t i = 0; i < length; ++i) {
        const float reference = routedInput[i + delay];
        const float reconstructed = separatedSum[i];
        const double delta = static_cast<double>(reference) - static_cast<double>(reconstructed);
        diffEnergy += delta * delta;
        referenceEnergy += static_cast<double>(reference) * static_cast<double>(reference);
    }

    ASSERT_GT(referenceEnergy, 0.0);
    const double normalizedError = std::sqrt(diffEnergy / referenceEnergy);
    std::cout << "HighFrequencyReconstruction normalized error: " << normalizedError << std::endl;
    EXPECT_LT(normalizedError, 0.02);
    EXPECT_NEAR(lastDetectedPitch, frequency, 3.0f);
}


