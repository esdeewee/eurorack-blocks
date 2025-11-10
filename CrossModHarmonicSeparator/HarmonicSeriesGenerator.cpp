#include "HarmonicSeriesGenerator.h"

#include <cmath>

namespace
{
constexpr float kNyquistFactor = 0.5f;
} // namespace

HarmonicSeriesGenerator::HarmonicSeriesGenerator(float sampleRate)
    : sample_rate_(sampleRate)
    , nyquist_(sampleRate * kNyquistFactor)
    , has_series_(false)
{
    reset();
}

void HarmonicSeriesGenerator::reset()
{
    fundamental_frequencies_.clear();
    prime_frequencies_.clear();
    composite_frequencies_.clear();
    prime_numbers_.clear();
    composite_numbers_.clear();
    has_series_ = false;
}

void HarmonicSeriesGenerator::update(float fundamentalHz)
{
    fundamental_frequencies_.clear();
    prime_frequencies_.clear();
    composite_frequencies_.clear();
    prime_numbers_.clear();
    composite_numbers_.clear();
    has_series_ = false;

    if (fundamentalHz <= 0.0f || fundamentalHz > nyquist_) {
        return;
    }

    std::size_t harmonicNumber = 1;

    while (true) {
        const float frequency = fundamentalHz * static_cast<float>(harmonicNumber);
        if (frequency > nyquist_) {
            break;
        }

        if (harmonicNumber == 1) {
            fundamental_frequencies_.push_back(frequency);
        } else if (isPrime(harmonicNumber)) {
            prime_numbers_.push_back(harmonicNumber);
            prime_frequencies_.push_back(frequency);
        } else {
            composite_numbers_.push_back(harmonicNumber);
            composite_frequencies_.push_back(frequency);
        }

        ++harmonicNumber;
    }

    has_series_ = !fundamental_frequencies_.empty();
}

bool HarmonicSeriesGenerator::hasSeries() const
{
    return has_series_;
}

const std::vector<float>& HarmonicSeriesGenerator::fundamentalFrequencies() const
{
    return fundamental_frequencies_;
}

const std::vector<float>& HarmonicSeriesGenerator::primeFrequencies() const
{
    return prime_frequencies_;
}

const std::vector<float>& HarmonicSeriesGenerator::compositeFrequencies() const
{
    return composite_frequencies_;
}

const std::vector<std::size_t>& HarmonicSeriesGenerator::primeNumbers() const
{
    return prime_numbers_;
}

const std::vector<std::size_t>& HarmonicSeriesGenerator::compositeNumbers() const
{
    return composite_numbers_;
}

bool HarmonicSeriesGenerator::isPrime(std::size_t n)
{
    if (n <= 1) {
        return false;
    }
    if (n == 2) {
        return true;
    }
    if ((n % 2) == 0) {
        return false;
    }

    const std::size_t limit = static_cast<std::size_t>(std::sqrt(static_cast<double>(n)));
    for (std::size_t i = 3; i <= limit; i += 2) {
        if ((n % i) == 0) {
            return false;
        }
    }
    return true;
}


