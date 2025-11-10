#include "HarmonicSpectralSeparator.h"

#include <algorithm>
#include <cmath>

namespace
{
constexpr float kPi = 3.14159265358979323846f;
} // namespace

HarmonicSpectralSeparator::HarmonicSpectralSeparator(float sampleRate,
                                                     std::size_t analysisSize,
                                                     std::size_t hopSize)
    : sample_rate_(sampleRate)
    , analysis_size_(analysisSize)
    , hop_size_(hopSize)
    , has_configuration_(false)
    , has_result_(false)
    , fundamental_hz_(0.0f)
    , total_energy_(0.0f)
    , fundamental_energy_(0.0f)
    , prime_energy_(0.0f)
    , composite_energy_(0.0f)
    , reconstruction_error_(0.0f)
{
    windowed_input_.resize(analysis_size_, 0.0f);
    spectrum_.resize(analysis_size_);
    fundamental_spectrum_.resize(analysis_size_);
    prime_spectrum_.resize(analysis_size_);
    composite_spectrum_.resize(analysis_size_);
    fundamental_time_domain_.resize(analysis_size_, 0.0f);
    prime_time_domain_.resize(analysis_size_, 0.0f);
    composite_time_domain_.resize(analysis_size_, 0.0f);
    ensureWindow();
    accumulator_.reserve(analysis_size_ * 2);
}

void HarmonicSpectralSeparator::reset()
{
    has_configuration_ = false;
    has_result_ = false;
    fundamental_hz_ = 0.0f;
    prime_numbers_.clear();
    composite_numbers_.clear();
    accumulator_.clear();
    std::fill(windowed_input_.begin(), windowed_input_.end(), 0.0f);
    clearSpectra();
    std::fill(fundamental_time_domain_.begin(), fundamental_time_domain_.end(), 0.0f);
    std::fill(prime_time_domain_.begin(), prime_time_domain_.end(), 0.0f);
    std::fill(composite_time_domain_.begin(), composite_time_domain_.end(), 0.0f);
    total_energy_ = fundamental_energy_ = prime_energy_ = composite_energy_ = 0.0f;
    reconstruction_error_ = 0.0f;
}

void HarmonicSpectralSeparator::setHarmonicData(float fundamentalHz,
                                                const std::vector<std::size_t>& primeNumbers,
                                                const std::vector<std::size_t>& compositeNumbers)
{
    if (fundamentalHz <= 0.0f || fundamentalHz * static_cast<float>(analysis_size_) / sample_rate_ >= static_cast<float>(analysis_size_ / 2)) {
        has_configuration_ = false;
        fundamental_hz_ = 0.0f;
        prime_numbers_.clear();
        composite_numbers_.clear();
        return;
    }

    fundamental_hz_ = fundamentalHz;
    prime_numbers_ = primeNumbers;
    composite_numbers_ = compositeNumbers;
    has_configuration_ = true;
}

bool HarmonicSpectralSeparator::processBuffer(const float* data, std::size_t length)
{
    if (data == nullptr || length == 0) {
        return false;
    }

    accumulator_.insert(accumulator_.end(), data, data + length);

    if (!has_configuration_ || accumulator_.size() < analysis_size_) {
        has_result_ = false;
        return false;
    }

    ensureWindow();

    const float* windowStart = accumulator_.data() + (accumulator_.size() - analysis_size_);
    for (std::size_t n = 0; n < analysis_size_; ++n) {
        windowed_input_[n] = windowStart[n] * window_[n];
    }

    total_energy_ = 0.0f;
    for (float sample : windowed_input_) {
        total_energy_ += sample * sample;
    }

    const float normFactor = -2.0f * kPi / static_cast<float>(analysis_size_);
    for (std::size_t k = 0; k < analysis_size_; ++k) {
        std::complex<float> sum {0.0f, 0.0f};
        for (std::size_t n = 0; n < analysis_size_; ++n) {
            const float angle = normFactor * static_cast<float>(k) * static_cast<float>(n);
            const float cosVal = std::cos(angle);
            const float sinVal = std::sin(angle);
            sum += std::complex<float> { windowed_input_[n] * cosVal, windowed_input_[n] * sinVal };
        }
        spectrum_[k] = sum;
    }

    clearSpectra();

    auto assignHarmonic = [&](std::size_t harmonicNumber, std::vector<std::complex<float>>& target) {
        if (harmonicNumber == 0) {
            return;
        }
        const float freq = fundamental_hz_ * static_cast<float>(harmonicNumber);
        const auto bin = static_cast<std::size_t>(std::llround(freq * static_cast<float>(analysis_size_) / sample_rate_));
        if (bin >= analysis_size_) {
            return;
        }
        assignBin(bin, target);
    };

    assignHarmonic(1, fundamental_spectrum_);
    for (std::size_t prime : prime_numbers_) {
        assignHarmonic(prime, prime_spectrum_);
    }
    for (std::size_t composite : composite_numbers_) {
        assignHarmonic(composite, composite_spectrum_);
    }

    computeInverse(fundamental_spectrum_, fundamental_time_domain_);
    computeInverse(prime_spectrum_, prime_time_domain_);
    computeInverse(composite_spectrum_, composite_time_domain_);

    auto computeEnergy = [](const std::vector<float>& buffer) {
        double sum = 0.0;
        for (float sample : buffer) {
            sum += static_cast<double>(sample) * static_cast<double>(sample);
        }
        return static_cast<float>(sum);
    };

    fundamental_energy_ = computeEnergy(fundamental_time_domain_);
    prime_energy_ = computeEnergy(prime_time_domain_);
    composite_energy_ = computeEnergy(composite_time_domain_);

    std::vector<float> reconstruction(analysis_size_, 0.0f);
    for (std::size_t n = 0; n < analysis_size_; ++n) {
        reconstruction[n] = fundamental_time_domain_[n]
                            + prime_time_domain_[n]
                            + composite_time_domain_[n];
    }

    double numerator = 0.0;
    double denominator = 0.0;
    for (std::size_t n = 0; n < analysis_size_; ++n) {
        const double diff = static_cast<double>(windowed_input_[n]) - static_cast<double>(reconstruction[n]);
        numerator += diff * diff;
        denominator += static_cast<double>(windowed_input_[n]) * static_cast<double>(windowed_input_[n]);
    }

    if (denominator <= 0.0) {
        reconstruction_error_ = 0.0f;
    } else {
        reconstruction_error_ = static_cast<float>(std::sqrt(numerator / denominator));
    }

    has_result_ = true;

    if (accumulator_.size() >= analysis_size_) {
        if (accumulator_.size() > analysis_size_) {
            accumulator_.erase(accumulator_.begin(), accumulator_.end() - analysis_size_ + hop_size_);
        } else {
            accumulator_.clear();
        }
    }

    return true;
}

bool HarmonicSpectralSeparator::hasResult() const
{
    return has_result_;
}

bool HarmonicSpectralSeparator::hasConfiguration() const
{
    return has_configuration_;
}

const std::vector<float>& HarmonicSpectralSeparator::windowedInput() const
{
    return windowed_input_;
}

const std::vector<float>& HarmonicSpectralSeparator::fundamentalTimeDomain() const
{
    return fundamental_time_domain_;
}

const std::vector<float>& HarmonicSpectralSeparator::primeTimeDomain() const
{
    return prime_time_domain_;
}

const std::vector<float>& HarmonicSpectralSeparator::compositeTimeDomain() const
{
    return composite_time_domain_;
}

float HarmonicSpectralSeparator::totalEnergy() const
{
    return total_energy_;
}

float HarmonicSpectralSeparator::fundamentalEnergy() const
{
    return fundamental_energy_;
}

float HarmonicSpectralSeparator::primeEnergy() const
{
    return prime_energy_;
}

float HarmonicSpectralSeparator::compositeEnergy() const
{
    return composite_energy_;
}

float HarmonicSpectralSeparator::reconstructionError() const
{
    return reconstruction_error_;
}

void HarmonicSpectralSeparator::ensureWindow()
{
    if (window_.size() == analysis_size_) {
        return;
    }
    window_.resize(analysis_size_);
    for (std::size_t n = 0; n < analysis_size_; ++n) {
        window_[n] = 0.5f * (1.0f - std::cos((2.0f * kPi * static_cast<float>(n)) / static_cast<float>(analysis_size_ - 1)));
    }
}

void HarmonicSpectralSeparator::clearSpectra()
{
    std::fill(fundamental_spectrum_.begin(), fundamental_spectrum_.end(), std::complex<float>{0.0f, 0.0f});
    std::fill(prime_spectrum_.begin(), prime_spectrum_.end(), std::complex<float>{0.0f, 0.0f});
    std::fill(composite_spectrum_.begin(), composite_spectrum_.end(), std::complex<float>{0.0f, 0.0f});
}

void HarmonicSpectralSeparator::assignBin(std::size_t bin, std::vector<std::complex<float>>& target)
{
    if (bin >= analysis_size_) {
        return;
    }

    target[bin] = spectrum_[bin];

    if (bin == 0 || bin == analysis_size_ / 2) {
        return;
    }

    const std::size_t mirror = analysis_size_ - bin;
    if (mirror < analysis_size_) {
        target[mirror] = spectrum_[mirror];
    }
}

void HarmonicSpectralSeparator::computeInverse(const std::vector<std::complex<float>>& spectrum,
                                               std::vector<float>& destination)
{
    const float norm = 1.0f / static_cast<float>(analysis_size_);
    const float angleFactor = 2.0f * kPi / static_cast<float>(analysis_size_);

    for (std::size_t n = 0; n < analysis_size_; ++n) {
        std::complex<float> sum {0.0f, 0.0f};
        for (std::size_t k = 0; k < analysis_size_; ++k) {
            const float angle = angleFactor * static_cast<float>(k) * static_cast<float>(n);
            const float cosVal = std::cos(angle);
            const float sinVal = std::sin(angle);
            const std::complex<float> basis { cosVal, sinVal };
            sum += spectrum[k] * basis;
        }
        destination[n] = norm * sum.real();
    }
}


