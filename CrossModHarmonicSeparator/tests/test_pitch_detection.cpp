#include <gtest/gtest.h>

#define CROSSMOD_TEST_ENV_HEADER "CrossModTestEnv.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "CrossModHarmonicSeparator/PitchDetector.h"
#include "helpers/TestSignalGenerators.h"

namespace
{
constexpr float kSampleRate = 48000.0f;
constexpr std::size_t kChunk = 48;

float amplitudeFromDb(float db)
{
    return std::pow(10.0f, db / 20.0f);
}

void feedSignal(PitchDetector& detector, const std::vector<float>& signal)
{
    std::size_t index = 0;
    while (index < signal.size()) {
        const std::size_t remaining = signal.size() - index;
        const std::size_t chunk = std::min<std::size_t>(kChunk, remaining);
        detector.processBuffer(signal.data() + index, chunk);
        index += chunk;
    }
}
} // namespace

TEST(PitchDetection, KnownFrequenciesAccuracy)
{
    PitchDetector detector(kSampleRate);

    const std::vector<float> frequencies { 100.0f, 200.0f, 440.0f, 1000.0f, 2000.0f };
    const std::vector<float> levelsDb { -40.0f, -20.0f, -10.0f, 0.0f };

    for (float freq : frequencies) {
        for (float db : levelsDb) {
            detector.reset();
            const float amplitude = amplitudeFromDb(db);
            auto signal = generateSine(freq, amplitude, kSampleRate, 2048);
            feedSignal(detector, signal);

            ASSERT_TRUE(detector.hasPitch()) << "Expected pitch detected for " << freq << " Hz at " << db << " dB";
            EXPECT_NEAR(detector.getCurrentPitchHz(), freq, 2.0f);
        }
    }
}

TEST(PitchDetection, HarmonicSignal)
{
    PitchDetector detector(kSampleRate);
    auto signal = generateSawtooth(200.0f, 0.5f, kSampleRate, 2048);
    feedSignal(detector, signal);

    ASSERT_TRUE(detector.hasPitch());
    EXPECT_NEAR(detector.getCurrentPitchHz(), 200.0f, 3.0f);
}

TEST(PitchDetection, SilenceHandling)
{
    PitchDetector detector(kSampleRate);

    auto tone = generateSine(440.0f, 0.5f, kSampleRate, 2048);
    feedSignal(detector, tone);
    ASSERT_TRUE(detector.hasPitch());
    EXPECT_NEAR(detector.getCurrentPitchHz(), 440.0f, 2.0f);

    auto silence = std::vector<float>(2048, 0.0f);
    feedSignal(detector, silence);

    EXPECT_FALSE(detector.hasPitch());
    EXPECT_NEAR(detector.getLastPitchHz(), 440.0f, 5.0f);
}

TEST(PitchDetection, RequiresSufficientSamples)
{
    PitchDetector detector(kSampleRate);

    auto shortTone = generateSine(440.0f, 0.5f, kSampleRate, 512);
    feedSignal(detector, shortTone);

    EXPECT_FALSE(detector.hasPitch());
    EXPECT_FLOAT_EQ(detector.getCurrentPitchHz(), 0.0f);
}

TEST(PitchDetection, SineSweepAccuracy)
{
    PitchDetector detector(kSampleRate);
    const float startHz = 80.0f;
    const float endHz = 2000.0f;

    auto sweep = generateSineSweep(startHz, endHz, 0.6f, kSampleRate, 48000);
    std::vector<float> errors;

    std::size_t index = 0;
    while (index < sweep.size()) {
        const std::size_t chunk = std::min<std::size_t>(kChunk, sweep.size() - index);
        detector.processBuffer(sweep.data() + index, chunk);

        if (detector.hasPitch()) {
            const float mid = static_cast<float>(index + chunk / 2);
            const float ratio = mid / static_cast<float>(sweep.size() - 1);
            const float expected = startHz + (endHz - startHz) * ratio;
            errors.push_back(std::fabs(detector.getCurrentPitchHz() - expected));
        }

        index += chunk;
    }

    ASSERT_FALSE(errors.empty());
    const float maxError = *std::max_element(errors.begin(), errors.end());
    EXPECT_LT(maxError, 25.0f);
}

TEST(PitchDetection, MultiToneRejection)
{
    PitchDetector detector(kSampleRate);
    std::vector<float> freqs { 440.0f, 660.0f };
    std::vector<float> amps { 0.4f, 0.4f };
    auto signal = generateMultiSine(freqs, amps, kSampleRate, 2048);

    feedSignal(detector, signal);

    ASSERT_TRUE(detector.hasPitch());
    EXPECT_NEAR(detector.getCurrentPitchHz(), 220.0f, 5.0f);
}

TEST(PitchDetection, LatencyWithinBudget)
{
    PitchDetector detector(kSampleRate);
    auto tone = generateSine(440.0f, 0.6f, kSampleRate, 4096);

    std::size_t index = 0;
    int buffersProcessed = 0;
    int firstDetection = -1;

    while (index < tone.size()) {
        const std::size_t chunk = std::min<std::size_t>(kChunk, tone.size() - index);
        detector.processBuffer(tone.data() + index, chunk);

        if (detector.hasPitch() && firstDetection == -1) {
            firstDetection = buffersProcessed;
        }

        index += chunk;
        ++buffersProcessed;
    }

    ASSERT_NE(firstDetection, -1);
    EXPECT_LE(firstDetection, 25);
}

