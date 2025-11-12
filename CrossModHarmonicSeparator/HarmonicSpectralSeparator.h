#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "../submodules/kissfft/kiss_fft.h"

class HarmonicSpectralSeparator
{
public:
    HarmonicSpectralSeparator(float sampleRate,
                              std::size_t analysisSize = 1024,
                              std::size_t hopSize = 48);
    ~HarmonicSpectralSeparator();

    HarmonicSpectralSeparator(const HarmonicSpectralSeparator&) = delete;
    HarmonicSpectralSeparator& operator=(const HarmonicSpectralSeparator&) = delete;
    HarmonicSpectralSeparator(HarmonicSpectralSeparator&&) = delete;
    HarmonicSpectralSeparator& operator=(HarmonicSpectralSeparator&&) = delete;

    void reset();
    void setHarmonicData(float fundamentalHz,
                         const std::vector<std::size_t>& primeNumbers,
                         const std::vector<std::size_t>& compositeNumbers);

    bool processBuffer(const float* data, std::size_t length);

    bool hasResult() const;
    bool hasConfiguration() const;

    const std::vector<float>& windowedInput() const;
    const std::vector<float>& fundamentalTimeDomain() const;
    const std::vector<float>& primeTimeDomain() const;
    const std::vector<float>& compositeTimeDomain() const;

    float totalEnergy() const;
    float fundamentalEnergy() const;
    float primeEnergy() const;
    float compositeEnergy() const;
    float reconstructionError() const;
    std::size_t analysisSize() const;
    std::size_t hopSize() const;
    const std::vector<float>& window() const;

private:
    void ensureWindow();
    void clearSpectra();
    void assignBin(std::size_t bin,
                   std::vector<std::complex<float>>& target,
                   std::vector<int>& ownership,
                   int category);
    void computeInverse(const std::vector<std::complex<float>>& spectrum,
                        std::vector<float>& destination);
    void allocateFft();
    void releaseFft();

    float sample_rate_;
    std::size_t analysis_size_;
    std::size_t hop_size_;

    bool has_configuration_;
    bool has_result_;

    float fundamental_hz_;
    std::vector<std::size_t> prime_numbers_;
    std::vector<std::size_t> composite_numbers_;

    std::vector<float> window_;
    std::vector<float> accumulator_;
    std::vector<float> windowed_input_;

    std::vector<std::complex<float>> spectrum_;
    std::vector<std::complex<float>> fundamental_spectrum_;
    std::vector<std::complex<float>> prime_spectrum_;
    std::vector<std::complex<float>> composite_spectrum_;

    std::vector<float> fundamental_time_domain_;
    std::vector<float> prime_time_domain_;
    std::vector<float> composite_time_domain_;

    kiss_fft_cfg forward_fft_;
    kiss_fft_cfg inverse_fft_;
    std::vector<kiss_fft_cpx> fft_input_;
    std::vector<kiss_fft_cpx> fft_output_;
    std::vector<kiss_fft_cpx> inverse_input_;
    std::vector<kiss_fft_cpx> inverse_output_;

    float total_energy_;
    float fundamental_energy_;
    float prime_energy_;
    float composite_energy_;
    float reconstruction_error_;
};


