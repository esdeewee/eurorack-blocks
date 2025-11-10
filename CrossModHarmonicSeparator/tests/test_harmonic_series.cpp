#include <gtest/gtest.h>

#define CROSSMOD_TEST_ENV_HEADER "CrossModTestEnv.h"

#include <algorithm>
#include <vector>

#include "HarmonicSeriesGenerator.h"

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


