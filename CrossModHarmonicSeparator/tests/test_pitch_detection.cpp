#include <gtest/gtest.h>

#define CROSSMOD_TEST_ENV_HEADER "CrossModTestEnv.h"

#include <cmath>
#include <vector>

#include "PitchDetector.h"
#include "helpers/TestSignalGenerators.h"

namespace
{
constexpr float kSampleRate = 48000.0f;
constexpr std::size_t kChunk = 48;

float amplitudeFromDb(float db)
{
    return std::pow(10.0f, db / 20.0f);
}

std::vector<float> generateSaw(float frequencyHz, float amplitude, float sampleRate, std::size_t sampleCount)
{
    std::vector<float> buffer(sampleCount);
    const float period = sampleRate / frequencyHz;

    for (std::size_t i = 0; i < sampleCount; ++i) {
        const float t = static_cast<float>(i);
        const float phase = std::fmod(t, period) / period;
        buffer[i] = amplitude * (2.0f * phase - 1.0f);
    }

    return buffer;
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
    auto signal = generateSaw(200.0f, 0.5f, kSampleRate, 2048);
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


