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



