#include <gtest/gtest.h>

#define CROSSMOD_TEST_ENV_HEADER "CrossModTestEnv.h"

#include <algorithm>
#include <random>
#include <vector>

#include "CrossModHarmonicSeparator/HarmonicSeriesGenerator.h"
#include "CrossModHarmonicSeparator/HarmonicSpectralSeparator.h"
#include "helpers/TestSignalGenerators.h"
#include "helpers/DSPTestUtils.h"

namespace
{
constexpr float kSampleRate = 48000.0f;
constexpr float kNyquist = kSampleRate * 0.5f;
} // namespace

TEST(HarmonicSeries, FrequencyCalculation)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(100.0f);

    ASSERT_TRUE(generator.hasSeries());

    const auto& fundamental = generator.fundamentalFrequencies();
    ASSERT_EQ(fundamental.size(), 1u);
    EXPECT_NEAR(fundamental.front(), 100.0f, 1e-4f);

    const auto& primeFreqs = generator.primeFrequencies();
    const auto& primeNums = generator.primeNumbers();
    ASSERT_GE(primeFreqs.size(), 4u);
    EXPECT_NEAR(primeFreqs[0], 200.0f, 1e-4f);
    EXPECT_NEAR(primeFreqs[1], 300.0f, 1e-4f);
    EXPECT_NEAR(primeFreqs[2], 500.0f, 1e-4f);
    EXPECT_NEAR(primeFreqs[3], 700.0f, 1e-4f);
    EXPECT_EQ(primeNums[0], 2u);
    EXPECT_EQ(primeNums[1], 3u);
    EXPECT_EQ(primeNums[2], 5u);
    EXPECT_EQ(primeNums[3], 7u);

    const auto& compositeFreqs = generator.compositeFrequencies();
    const auto& compositeNums = generator.compositeNumbers();
    ASSERT_GE(compositeFreqs.size(), 3u);
    EXPECT_NEAR(compositeFreqs[0], 400.0f, 1e-4f);
    EXPECT_NEAR(compositeFreqs[1], 600.0f, 1e-4f);
    EXPECT_NEAR(compositeFreqs[2], 800.0f, 1e-4f);
    EXPECT_EQ(compositeNums[0], 4u);
    EXPECT_EQ(compositeNums[1], 6u);
    EXPECT_EQ(compositeNums[2], 8u);
}

TEST(HarmonicSeries, PrimeCategorization)
{
    std::vector<std::size_t> expectedPrimes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47
    };

    for (std::size_t n = 1; n <= 50; ++n) {
        const bool isPrime = HarmonicSeriesGenerator::isPrime(n);
        const bool inList = std::find(expectedPrimes.begin(), expectedPrimes.end(), n) != expectedPrimes.end();

        if (n == 1) {
            EXPECT_FALSE(isPrime);
        } else if (inList) {
            EXPECT_TRUE(isPrime) << "Expected prime: " << n;
        } else {
            EXPECT_FALSE(isPrime) << "Expected composite: " << n;
        }
    }
}

TEST(HarmonicSeries, NyquistLimitHandling)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(10000.0f);

    ASSERT_TRUE(generator.hasSeries());
    const auto& fundamental = generator.fundamentalFrequencies();
    ASSERT_EQ(fundamental.size(), 1u);
    EXPECT_NEAR(fundamental.front(), 10000.0f, 1e-4f);

    const auto& primes = generator.primeFrequencies();
    ASSERT_EQ(primes.size(), 1u);
    EXPECT_NEAR(primes.front(), 20000.0f, 1e-4f);
    EXPECT_LT(primes.front(), kNyquist + 1.0f);
    EXPECT_TRUE(generator.compositeFrequencies().empty());

    generator.update(25000.0f);
    EXPECT_FALSE(generator.hasSeries());
    EXPECT_TRUE(generator.fundamentalFrequencies().empty());
}

TEST(HarmonicSeries, InvalidFundamental)
{
    HarmonicSeriesGenerator generator(kSampleRate);

    generator.update(0.0f);
    EXPECT_FALSE(generator.hasSeries());

    generator.update(-100.0f);
    EXPECT_FALSE(generator.hasSeries());
}

TEST(HarmonicSeries, LargeHarmonicCount)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(55.0f);

    ASSERT_TRUE(generator.hasSeries());

    std::size_t expectedPrimeCount = 0;
    std::size_t expectedCompositeCount = 0;
    std::size_t harmonic = 2;
    while (true) {
        const float freq = 55.0f * static_cast<float>(harmonic);
        if (freq > kNyquist) {
            break;
        }
        if (HarmonicSeriesGenerator::isPrime(harmonic)) {
            ++expectedPrimeCount;
        } else {
            ++expectedCompositeCount;
        }
        ++harmonic;
    }

    EXPECT_EQ(generator.primeNumbers().size(), expectedPrimeCount);
    EXPECT_EQ(generator.compositeNumbers().size(), expectedCompositeCount);
}

TEST(HarmonicSeries, JitterRobustness)
{
    HarmonicSeriesGenerator generator(kSampleRate);

    std::vector<float> fundamentals { 100.0f, 101.0f, 99.5f, 100.5f, 100.2f };
    for (float f0 : fundamentals) {
        generator.update(f0);
        ASSERT_TRUE(generator.hasSeries());
        ASSERT_FALSE(generator.fundamentalFrequencies().empty());
        EXPECT_NEAR(generator.fundamentalFrequencies().front(), f0, 0.5f);
    }
}

TEST(HarmonicSeries, RatiosAlignWithIntegers)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(120.0f);
    ASSERT_TRUE(generator.hasSeries());

    const float fundamental = generator.fundamentalFrequencies().front();

    auto validate = [&](const std::vector<float>& freqs, const std::vector<std::size_t>& numbers) {
        ASSERT_EQ(freqs.size(), numbers.size());
        for (std::size_t i = 0; i < freqs.size(); ++i) {
            const float ratio = freqs[i] / fundamental;
            EXPECT_NEAR(ratio, static_cast<float>(numbers[i]), 0.05f);
        }
    };

    validate(generator.primeFrequencies(), generator.primeNumbers());
    validate(generator.compositeFrequencies(), generator.compositeNumbers());
}

TEST(HarmonicSpectralSeparator, FundamentalCapturesEnergy)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(100.0f);

    HarmonicSpectralSeparator separator(kSampleRate);
    separator.setHarmonicData(100.0f, generator.primeNumbers(), generator.compositeNumbers());

    auto sine = generateSine(100.0f, 0.7f, kSampleRate, 1024);
    ASSERT_TRUE(separator.processBuffer(sine.data(), sine.size()));
    ASSERT_TRUE(separator.hasResult());

    const float total = separator.totalEnergy();
    ASSERT_GT(total, 0.0f);
    EXPECT_GT(separator.fundamentalEnergy() / total, 0.6f);
    EXPECT_LT(separator.primeEnergy() / total, 0.2f);
    EXPECT_LT(separator.compositeEnergy() / total, 0.2f);
}

TEST(HarmonicSpectralSeparator, PrimeCompositeIsolation)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(100.0f);

    HarmonicSpectralSeparator separator(kSampleRate);
    separator.setHarmonicData(100.0f, generator.primeNumbers(), generator.compositeNumbers());

    std::vector<std::size_t> primeHarmonics { 2, 3, 5 };
    auto primeSignal = generateHarmonicSum(100.0f, primeHarmonics, 0.6f, kSampleRate, 1024);

    ASSERT_TRUE(separator.processBuffer(primeSignal.data(), primeSignal.size()));
    const float totalPrime = separator.totalEnergy();
    EXPECT_GT(separator.primeEnergy() / totalPrime, 0.5f);
    EXPECT_LT(separator.compositeEnergy() / totalPrime, 0.35f);

    HarmonicSpectralSeparator compositeSeparator(kSampleRate);
    compositeSeparator.setHarmonicData(100.0f, generator.primeNumbers(), generator.compositeNumbers());

    std::vector<std::size_t> compositeHarmonics { 4, 6, 8 };
    auto compositeSignal = generateHarmonicSum(100.0f, compositeHarmonics, 0.6f, kSampleRate, 1024);
    ASSERT_TRUE(compositeSeparator.processBuffer(compositeSignal.data(), compositeSignal.size()));
    const float totalComposite = compositeSeparator.totalEnergy();
    EXPECT_GT(compositeSeparator.compositeEnergy() / totalComposite, 0.5f);
    EXPECT_LT(compositeSeparator.primeEnergy() / totalComposite, 0.35f);
}

TEST(HarmonicSpectralSeparator, ReconstructionAccuracy)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(120.0f);

    HarmonicSpectralSeparator separator(kSampleRate);
    separator.setHarmonicData(120.0f, generator.primeNumbers(), generator.compositeNumbers());

    std::vector<std::size_t> harmonics { 1, 2, 3, 4, 5 };
    auto complexSignal = generateHarmonicSum(120.0f, harmonics, 0.7f, kSampleRate, 1024);
    ASSERT_TRUE(separator.processBuffer(complexSignal.data(), complexSignal.size()));
    ASSERT_TRUE(separator.hasResult());

    EXPECT_LT(separator.reconstructionError(), 0.05f);
}

TEST(HarmonicSpectralSeparator, ReconstructsLowFundamental)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(60.0f);
    ASSERT_TRUE(generator.hasSeries());

    HarmonicSpectralSeparator separator(kSampleRate);
    separator.setHarmonicData(60.0f, generator.primeNumbers(), generator.compositeNumbers());

    std::vector<std::size_t> harmonics { 1 };
    harmonics.insert(harmonics.end(), generator.primeNumbers().begin(), generator.primeNumbers().end());
    harmonics.insert(harmonics.end(), generator.compositeNumbers().begin(), generator.compositeNumbers().end());

    auto signal = generateHarmonicSum(60.0f, harmonics, 0.6f, kSampleRate, 4096);
    ASSERT_TRUE(separator.processBuffer(signal.data(), signal.size()));
    ASSERT_TRUE(separator.hasResult());

    EXPECT_LT(separator.reconstructionError(), 0.05f);
}

TEST(HarmonicSpectralSeparator, ReconstructsHighFundamental)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(4000.0f);
    ASSERT_TRUE(generator.hasSeries());

    HarmonicSpectralSeparator separator(kSampleRate);
    separator.setHarmonicData(4000.0f, generator.primeNumbers(), generator.compositeNumbers());

    std::vector<std::size_t> harmonics { 1 };
    harmonics.insert(harmonics.end(), generator.primeNumbers().begin(), generator.primeNumbers().end());
    harmonics.insert(harmonics.end(), generator.compositeNumbers().begin(), generator.compositeNumbers().end());

    auto signal = generateHarmonicSum(4000.0f, harmonics, 0.6f, kSampleRate, 4096);
    ASSERT_TRUE(separator.processBuffer(signal.data(), signal.size()));
    ASSERT_TRUE(separator.hasResult());

    EXPECT_LT(separator.reconstructionError(), 0.05f);
}

TEST(HarmonicSpectralSeparator, NoiseInjectionTolerance)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(150.0f);

    HarmonicSpectralSeparator baseSeparator(kSampleRate);
    baseSeparator.setHarmonicData(150.0f, generator.primeNumbers(), generator.compositeNumbers());
    auto baseSignal = generateHarmonicSum(150.0f, generator.primeNumbers(), 0.6f, kSampleRate, 1024);
    ASSERT_TRUE(baseSeparator.processBuffer(baseSignal.data(), baseSignal.size()));
    const float baseRatio = baseSeparator.primeEnergy() / baseSeparator.totalEnergy();

    HarmonicSpectralSeparator noisySeparator(kSampleRate);
    noisySeparator.setHarmonicData(150.0f, generator.primeNumbers(), generator.compositeNumbers());

    std::vector<float> noisySignal = baseSignal;
    std::mt19937 rng(42);
    std::normal_distribution<float> dist(0.0f, 0.001f);
    for (float& sample : noisySignal) {
        sample += dist(rng);
    }

    ASSERT_TRUE(noisySeparator.processBuffer(noisySignal.data(), noisySignal.size()));
    const float noisyRatio = noisySeparator.primeEnergy() / noisySeparator.totalEnergy();

    EXPECT_NEAR(noisyRatio, baseRatio, 0.02f);
}

TEST(HarmonicSpectralSeparator, NyquistExclusion)
{
    HarmonicSeriesGenerator generator(kSampleRate);
    generator.update(10000.0f);

    HarmonicSpectralSeparator separator(kSampleRate);
    separator.setHarmonicData(10000.0f, generator.primeNumbers(), generator.compositeNumbers());

    auto signal = generateSine(10000.0f, 0.6f, kSampleRate, 1024);
    ASSERT_TRUE(separator.processBuffer(signal.data(), signal.size()));
    ASSERT_TRUE(separator.hasResult());

    EXPECT_GT(separator.fundamentalEnergy(), 0.0f);
    EXPECT_FLOAT_EQ(separator.compositeEnergy(), 0.0f);
}

TEST(HarmonicSpectralSeparator, ResetWhenNoConfiguration)
{
    HarmonicSpectralSeparator separator(kSampleRate);
    auto signal = generateSine(440.0f, 0.6f, kSampleRate, 1024);

    EXPECT_FALSE(separator.processBuffer(signal.data(), signal.size()));
    EXPECT_FALSE(separator.hasResult());
}



