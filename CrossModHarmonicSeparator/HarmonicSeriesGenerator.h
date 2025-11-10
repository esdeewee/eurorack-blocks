#pragma once

#include <cstddef>
#include <vector>

class HarmonicSeriesGenerator
{
public:
    explicit HarmonicSeriesGenerator(float sampleRate);

    void reset();
    void update(float fundamentalHz);

    bool hasSeries() const;

    const std::vector<float>& fundamentalFrequencies() const;
    const std::vector<float>& primeFrequencies() const;
    const std::vector<float>& compositeFrequencies() const;

    const std::vector<std::size_t>& primeNumbers() const;
    const std::vector<std::size_t>& compositeNumbers() const;

    static bool isPrime(std::size_t n);

private:
    float sample_rate_;
    float nyquist_;
    bool has_series_;

    std::vector<float> fundamental_frequencies_;
    std::vector<float> prime_frequencies_;
    std::vector<float> composite_frequencies_;

    std::vector<std::size_t> prime_numbers_;
    std::vector<std::size_t> composite_numbers_;
};


