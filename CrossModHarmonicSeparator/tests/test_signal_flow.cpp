#include <gtest/gtest.h>

#define CROSSMOD_TEST_ENV_HEADER "CrossModTestEnv.h"
#define CROSSMOD_DEBUG_ENABLED 1

#include <array>
#include <vector>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

#include "CrossModHarmonicSeparator.h"
#include "helpers/DSPTestUtils.h"
#include "helpers/TestSignalGenerators.h"

namespace
{

constexpr float kSampleRate = 48000.0f;

// Helper function to analyze energy preservation and spectral separation for different frequencies
void analyzeEnergyPreservationDetailed(const char* scenarioName,
                                      float frequency,
                                      float amplitude,
                                      std::size_t durationSamples = 48000)
{
    CrossModHarmonicSeparator module;
    module.init();

    auto source = generateSine(frequency, amplitude, kSampleRate, durationSamples);

    std::vector<float> routedInput;
    std::vector<float> separatedSum;
    std::vector<float> fundamentalOnly;
    std::vector<float> primeOnly;
    std::vector<float> compositeOnly;

    double detectedPitchSum = 0.0;
    int detectedPitchCount = 0;
    bool hasSpectralSeparation = false;

    std::size_t cursor = 0;
    while (cursor < source.size()) {
        for (std::size_t i = 0; i < erb_BUFFER_SIZE; ++i) {
            const float sample = (cursor < source.size()) ? source[cursor++] : 0.0f;
            module.ui.audio_in[i] = sample;
            routedInput.push_back(sample);
        }

        module.process();
        
        if (module.has_detected_pitch) {
            detectedPitchSum += module.detected_pitch_hz;
            detectedPitchCount++;
        }
        if (module.has_spectral_separation) {
            hasSpectralSeparation = true;
        }

        for (std::size_t i = 0; i < erb_BUFFER_SIZE; ++i) {
            const float fund = module.ui.fundamental_debug[i];
            const float prime = module.ui.prime_debug[i];
            const float comp = module.ui.composite_debug[i];
            const float separated = fund + prime + comp;
            
            separatedSum.push_back(separated);
            fundamentalOnly.push_back(fund);
            primeOnly.push_back(prime);
            compositeOnly.push_back(comp);
        }
    }

    const std::size_t analysis = module.spectral_separator.analysisSize();
    const std::size_t hop = module.spectral_separator.hopSize();
    const std::size_t delay = (analysis > hop) ? (analysis - hop) : 0;

    ASSERT_GT(routedInput.size(), delay + 1024) << scenarioName << ": Insufficient input samples";
    ASSERT_GT(separatedSum.size(), delay + 1024) << scenarioName << ": Insufficient output samples";

    const std::size_t length = std::min(routedInput.size() - delay, separatedSum.size() - delay);
    
    // Calculate energies
    double energyInput = 0.0;
    double energySum = 0.0;
    double energyFundamental = 0.0;
    double energyPrime = 0.0;
    double energyComposite = 0.0;
    
    for (std::size_t i = 0; i < length; ++i) {
        const std::size_t inputIdx = delay + i;
        const std::size_t outputIdx = delay + i;
        if (inputIdx < routedInput.size() && outputIdx < separatedSum.size()) {
            const double inVal = static_cast<double>(routedInput[inputIdx]);
            const double sumVal = static_cast<double>(separatedSum[outputIdx]);
            const double fundVal = static_cast<double>(fundamentalOnly[outputIdx]);
            const double primeVal = static_cast<double>(primeOnly[outputIdx]);
            const double compVal = static_cast<double>(compositeOnly[outputIdx]);
            
            energyInput += inVal * inVal;
            energySum += sumVal * sumVal;
            energyFundamental += fundVal * fundVal;
            energyPrime += primeVal * primeVal;
            energyComposite += compVal * compVal;
        }
    }

    const double energyRatio = (energyInput > 0.0) ? (energySum / energyInput) : 0.0;
    const double avgDetectedPitch = (detectedPitchCount > 0) ? (detectedPitchSum / detectedPitchCount) : 0.0;
    const double pitchError = std::abs(avgDetectedPitch - frequency);
    const double pitchErrorPercent = (frequency > 0.0) ? (pitchError / frequency * 100.0) : 0.0;
    
    // Calculate energy distribution
    const double totalSeparatedEnergy = energyFundamental + energyPrime + energyComposite;
    const double fundPercent = (totalSeparatedEnergy > 0.0) ? (energyFundamental / totalSeparatedEnergy * 100.0) : 0.0;
    const double primePercent = (totalSeparatedEnergy > 0.0) ? (energyPrime / totalSeparatedEnergy * 100.0) : 0.0;
    const double compPercent = (totalSeparatedEnergy > 0.0) ? (energyComposite / totalSeparatedEnergy * 100.0) : 0.0;
    
    std::cout << "\n=== [" << scenarioName << "] Detailed Analysis ===\n";
    std::cout << "Frequency: " << frequency << " Hz, Amplitude: " << amplitude << "\n";
    std::cout << "Detected pitch: " << avgDetectedPitch << " Hz (error: " << pitchError << " Hz, " << pitchErrorPercent << "%)\n";
    std::cout << "Pitch detection count: " << detectedPitchCount << "\n";
    std::cout << "Has spectral separation: " << (hasSpectralSeparation ? "yes" : "no") << "\n";
    std::cout << "\nEnergy Analysis:\n";
    std::cout << "  Input energy: " << energyInput << "\n";
    std::cout << "  Output energy (sum): " << energySum << "\n";
    std::cout << "  Energy ratio: " << energyRatio << "\n";
    std::cout << "  Energy loss/gain: " << ((energyInput - energySum) / energyInput * 100.0) << "%\n";
    std::cout << "\nSpectral Separation Distribution:\n";
    std::cout << "  Fundamental: " << energyFundamental << " (" << fundPercent << "%)\n";
    std::cout << "  Prime: " << energyPrime << " (" << primePercent << "%)\n";
    std::cout << "  Composite: " << energyComposite << " (" << compPercent << "%)\n";
    std::cout << "  Total separated: " << totalSeparatedEnergy << "\n";
    std::cout << "========================================\n";
}

// Helper function to test energy preservation for different scenarios
void testEnergyPreservationScenario(const char* scenarioName,
                                    float frequency,
                                    float amplitude,
                                    std::size_t durationSamples = 48000)
{
    CrossModHarmonicSeparator module;
    module.init();

    auto source = generateSine(frequency, amplitude, kSampleRate, durationSamples);

    std::vector<float> routedInput;
    std::vector<float> separatedSum;

    std::size_t cursor = 0;
    while (cursor < source.size()) {
        for (std::size_t i = 0; i < erb_BUFFER_SIZE; ++i) {
            const float sample = (cursor < source.size()) ? source[cursor++] : 0.0f;
            module.ui.audio_in[i] = sample;
            routedInput.push_back(sample);
        }

        module.process();

        for (std::size_t i = 0; i < erb_BUFFER_SIZE; ++i) {
            const float separated = module.ui.fundamental_debug[i]
                                    + module.ui.prime_debug[i]
                                    + module.ui.composite_debug[i];
            separatedSum.push_back(separated);
        }
    }

    ASSERT_TRUE(module.has_spectral_separation) << scenarioName << ": Spectral separation should be active";

    const std::size_t analysis = module.spectral_separator.analysisSize();
    const std::size_t hop = module.spectral_separator.hopSize();
    const std::size_t delay = (analysis > hop) ? (analysis - hop) : 0;

    ASSERT_GT(routedInput.size(), delay + 1024) << scenarioName << ": Insufficient input samples";
    ASSERT_GT(separatedSum.size(), delay + 1024) << scenarioName << ": Insufficient output samples";

    const std::size_t length = std::min(routedInput.size() - delay, separatedSum.size() - delay);
    double energyFundamental = 0.0;
    double energySum = 0.0;
    for (std::size_t i = 0; i < length; ++i) {
        const std::size_t inputIdx = delay + i;
        const std::size_t outputIdx = delay + i;
        if (inputIdx < routedInput.size() && outputIdx < separatedSum.size()) {
            energyFundamental += static_cast<double>(routedInput[inputIdx]) * static_cast<double>(routedInput[inputIdx]);
            energySum += static_cast<double>(separatedSum[outputIdx]) * static_cast<double>(separatedSum[outputIdx]);
        }
    }

    const double energyRatio = (energyFundamental > 0.0) ? (energySum / energyFundamental) : 0.0;
    const double energyLossPercent = (energyFundamental > 0.0) ? ((energyFundamental - energySum) / energyFundamental * 100.0) : 0.0;
    
    std::cout << "\n[" << scenarioName << "] Frequency: " << frequency << " Hz, Amplitude: " << amplitude
              << ", Energy ratio: " << energyRatio << ", Loss: " << energyLossPercent << "%\n";

    // Verify energy preservation (should be within 2% with compensation factor 0.92)
    EXPECT_NEAR(energySum, energyFundamental, energyFundamental * 0.02)
        << scenarioName << ": Energy mismatch - Input: " << energyFundamental 
        << ", Output: " << energySum << ", Ratio: " << energyRatio;
    
    // Verify that the energy ratio is close to 1.0 (within 2%)
    EXPECT_NEAR(energyRatio, 1.0, 0.02)
        << scenarioName << ": Energy ratio should be close to 1.0 with compensation factor";
}

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

    const std::size_t length = std::min(routedInput.size() - delay, separatedSum.size() - delay);
    double energyFundamental = 0.0;
    double energySum = 0.0;
    for (std::size_t i = 0; i < length; ++i) {
        const std::size_t inputIdx = delay + i;
        const std::size_t outputIdx = delay + i;
        if (inputIdx < routedInput.size() && outputIdx < separatedSum.size()) {
            energyFundamental += static_cast<double>(routedInput[inputIdx]) * static_cast<double>(routedInput[inputIdx]);
            energySum += static_cast<double>(separatedSum[outputIdx]) * static_cast<double>(separatedSum[outputIdx]);
        }
    }

    EXPECT_NEAR(energySum, energyFundamental, energyFundamental * 0.02);
    
    // Debug output if test fails
    if (std::abs(energySum - energyFundamental) > energyFundamental * 0.02) {
        std::cout << "\n=== Energy Debug Info ===\n";
        std::cout << "Input energy: " << energyFundamental << "\n";
        std::cout << "Output energy: " << energySum << "\n";
        std::cout << "Energy ratio: " << (energySum / energyFundamental) << "\n";
        std::cout << "Energy loss: " << ((energyFundamental - energySum) / energyFundamental * 100.0) << "%\n";
        
        if (!module.normalization_pattern.empty()) {
            std::cout << "\nNormalization pattern stats:\n";
            double minNorm = *std::min_element(module.normalization_pattern.begin(), module.normalization_pattern.end());
            double maxNorm = *std::max_element(module.normalization_pattern.begin(), module.normalization_pattern.end());
            double avgNorm = std::accumulate(module.normalization_pattern.begin(), module.normalization_pattern.end(), 0.0) / module.normalization_pattern.size();
            std::cout << "  Min: " << minNorm << ", Max: " << maxNorm << ", Avg: " << avgNorm << "\n";
            std::cout << "  Pattern size: " << module.normalization_pattern.size() << "\n";
            
            #ifdef CROSSMOD_DEBUG_ENABLED
            std::cout << "\nCROSSMOD_DEBUG_ENABLED is defined\n";
            if (!module.debug_window_sum_pattern.empty()) {
                std::cout << "\nCOLA window sum pattern stats:\n";
                double minSum = *std::min_element(module.debug_window_sum_pattern.begin(), module.debug_window_sum_pattern.end());
                double maxSum = *std::max_element(module.debug_window_sum_pattern.begin(), module.debug_window_sum_pattern.end());
                double avgSum = std::accumulate(module.debug_window_sum_pattern.begin(), module.debug_window_sum_pattern.end(), 0.0) / module.debug_window_sum_pattern.size();
                std::cout << "  Min: " << minSum << ", Max: " << maxSum << ", Avg: " << avgSum << "\n";
                std::cout << "  COLA variation: " << ((maxSum - minSum) / avgSum * 100.0) << "%\n";
                std::cout << "  Pattern size: " << module.debug_window_sum_pattern.size() << "\n";
                
                // Calculate expected normalization if we divide by sum(w) instead of sum(w^2)
                if (avgSum > 0.0 && avgNorm > 0.0) {
                    double expected_energy_ratio = (avgSum * avgSum) / avgNorm;
                    std::cout << "  Expected energy ratio if dividing by sum(w): " << expected_energy_ratio << "\n";
                }
            } else {
                std::cout << "\nCOLA window sum pattern: EMPTY (not computed)\n";
            }
            
            if (!module.debug_overlap_count_pattern.empty()) {
                int minOverlap = *std::min_element(module.debug_overlap_count_pattern.begin(), module.debug_overlap_count_pattern.end());
                int maxOverlap = *std::max_element(module.debug_overlap_count_pattern.begin(), module.debug_overlap_count_pattern.end());
                double avgOverlap = std::accumulate(module.debug_overlap_count_pattern.begin(), module.debug_overlap_count_pattern.end(), 0.0) / module.debug_overlap_count_pattern.size();
                std::cout << "\nOverlap count pattern stats:\n";
                std::cout << "  Min: " << minOverlap << ", Max: " << maxOverlap << ", Avg: " << avgOverlap << "\n";
                std::cout << "  Pattern size: " << module.debug_overlap_count_pattern.size() << "\n";
            } else {
                std::cout << "\nOverlap count pattern: EMPTY (not computed)\n";
            }
            
            // Calculate theoretical compensation factor based on window sum and normalization pattern
            if (!module.debug_window_sum_pattern.empty() && !module.normalization_pattern.empty()) {
                // For energy preservation in overlap-add with IFFT normalization:
                // - IFFT outputs: x[n] / N
                // - After windowing: (x[n] / N) * w[n]
                // - After overlap-add: sum((x[k] / N) * w[k]) = (1/N) * sum(x[k] * w[k])
                // - Energy: E_overlap = (1/N^2) * sum((sum(x[k] * w[k]))^2)
                // - For a constant signal x: E_overlap = (x/N)^2 * (sum(w[k]))^2
                // - Dividing by sum(w^2): E_output = (x/N)^2 * (sum(w[k]))^2 / sum(w^2)
                // - Original energy: E_original = x^2
                // - So: E_output / E_original = (1/N^2) * (sum(w[k]))^2 / sum(w^2)
                // - Compensation factor = sqrt(E_original / E_output) = N * sqrt(sum(w^2)) / sum(w[k])
                double avgWindowSum = std::accumulate(module.debug_window_sum_pattern.begin(), module.debug_window_sum_pattern.end(), 0.0) / module.debug_window_sum_pattern.size();
                double avgNorm = std::accumulate(module.normalization_pattern.begin(), module.normalization_pattern.end(), 0.0) / module.normalization_pattern.size();
                double N = static_cast<double>(module.spectral_separator.analysisSize());
                double theoretical_compensation = N * std::sqrt(avgNorm) / avgWindowSum;
                std::cout << "\nTheoretical compensation factor analysis:\n";
                std::cout << "  N (analysis size): " << N << "\n";
                std::cout << "  Avg window sum: " << avgWindowSum << "\n";
                std::cout << "  Avg normalization (sum w^2): " << avgNorm << "\n";
                std::cout << "  Theoretical compensation: " << theoretical_compensation << "\n";
                std::cout << "  Empirical compensation (from energy ratio): " << (1.0 / std::sqrt(energySum / energyFundamental)) << "\n";
            }
            
            // Show energy measurements from debug members
            std::cout << "\nDebug energy measurements:\n";
            std::cout << "  Pre-normalization energy: " << module.debug_pre_normalization_energy << "\n";
            std::cout << "  Post-normalization energy: " << module.debug_post_normalization_energy << "\n";
            std::cout << "  Input energy (current buffer): " << module.debug_input_energy << "\n";
            if (module.debug_pre_normalization_energy > 0.0) {
                std::cout << "  Normalization energy ratio: " << (module.debug_post_normalization_energy / module.debug_pre_normalization_energy) << "\n";
            }
            if (!module.debug_normalization_values.empty()) {
                double minNormVal = *std::min_element(module.debug_normalization_values.begin(), module.debug_normalization_values.end());
                double maxNormVal = *std::max_element(module.debug_normalization_values.begin(), module.debug_normalization_values.end());
                double avgNormVal = std::accumulate(module.debug_normalization_values.begin(), module.debug_normalization_values.end(), 0.0) / module.debug_normalization_values.size();
                std::cout << "  Normalization values used: Min=" << minNormVal << ", Max=" << maxNormVal << ", Avg=" << avgNormVal << "\n";
            }
            #else
            std::cout << "\nDebug patterns not available (CROSSMOD_DEBUG_ENABLED not defined)\n";
            #endif
        }
        
        std::cout << "\nAnalysis size: " << module.spectral_separator.analysisSize() << "\n";
        std::cout << "Hop size: " << module.spectral_separator.hopSize() << "\n";
        std::cout << "Overlap ratio: " << (1.0 - static_cast<double>(module.spectral_separator.hopSize()) / module.spectral_separator.analysisSize()) * 100.0 << "%\n";
        std::cout << "Has spectral separation: " << (module.has_spectral_separation ? "yes" : "no") << "\n";
        std::cout << "Fundamental overlap size: " << module.fundamental_overlap.size() << "\n";
        std::cout << "========================\n";
    }
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

    const std::size_t length = std::min(routedInput.size() - delay, separatedSum.size() - delay);
    double diffEnergy = 0.0;
    double referenceEnergy = 0.0;
    for (std::size_t i = 0; i < length; ++i) {
        const std::size_t inputIdx = delay + i;
        const std::size_t outputIdx = delay + i;
        if (inputIdx < routedInput.size() && outputIdx < separatedSum.size()) {
            const float reference = routedInput[inputIdx];
            const float reconstructed = separatedSum[outputIdx];
            const double delta = static_cast<double>(reference) - static_cast<double>(reconstructed);
            diffEnergy += delta * delta;
            referenceEnergy += static_cast<double>(reference) * static_cast<double>(reference);
        }
    }

    ASSERT_GT(referenceEnergy, 0.0);
    const double normalizedError = std::sqrt(diffEnergy / referenceEnergy);
    std::cout << "HighFrequencyReconstruction normalized error: " << normalizedError << std::endl;
    EXPECT_LT(normalizedError, 0.02);
    EXPECT_NEAR(lastDetectedPitch, frequency, 3.0f);
}

// Detailed analysis tests to understand frequency dependency
TEST(SignalFlow, DetailedAnalysisLowFrequency)
{
    analyzeEnergyPreservationDetailed("LowFreq_Analysis", 60.0f, 0.5f, 24000);
}

TEST(SignalFlow, DetailedAnalysisMidFrequency)
{
    analyzeEnergyPreservationDetailed("MidFreq_Analysis", 440.0f, 0.5f, 24000);
}

TEST(SignalFlow, DetailedAnalysisHighFrequency)
{
    analyzeEnergyPreservationDetailed("HighFreq_Analysis", 2000.0f, 0.5f, 24000);
}

TEST(SignalFlow, DetailedAnalysisVeryLowFrequency)
{
    analyzeEnergyPreservationDetailed("VeryLowFreq_Analysis", 40.0f, 0.5f, 24000);
}

TEST(SignalFlow, EnergyPreservationLowFrequency)
{
    testEnergyPreservationScenario("LowFrequency", 60.0f, 0.5f);
}

TEST(SignalFlow, EnergyPreservationMidFrequency)
{
    testEnergyPreservationScenario("MidFrequency", 440.0f, 0.5f);
}

TEST(SignalFlow, EnergyPreservationHighFrequency)
{
    testEnergyPreservationScenario("HighFrequency", 2000.0f, 0.5f);
}

TEST(SignalFlow, EnergyPreservationVeryLowFrequency)
{
    testEnergyPreservationScenario("VeryLowFrequency", 40.0f, 0.5f);
}

TEST(SignalFlow, EnergyPreservationVeryHighFrequency)
{
    testEnergyPreservationScenario("VeryHighFrequency", 5000.0f, 0.4f);
}

TEST(SignalFlow, EnergyPreservationLowAmplitude)
{
    testEnergyPreservationScenario("LowAmplitude", 440.0f, 0.1f);
}

TEST(SignalFlow, EnergyPreservationHighAmplitude)
{
    testEnergyPreservationScenario("HighAmplitude", 440.0f, 0.9f);
}

TEST(SignalFlow, EnergyPreservationMediumAmplitude)
{
    testEnergyPreservationScenario("MediumAmplitude", 440.0f, 0.5f);
}

TEST(SignalFlow, EnergyPreservationMultipleFrequencies)
{
    // Test with multiple frequencies to verify consistency
    std::vector<float> frequencies { 100.0f, 440.0f, 1000.0f, 2000.0f };
    std::vector<float> amplitudes { 0.3f, 0.3f, 0.2f, 0.2f };
    
    for (std::size_t i = 0; i < frequencies.size(); ++i) {
        std::string scenarioName = "MultiFreq_" + std::to_string(static_cast<int>(frequencies[i]));
        testEnergyPreservationScenario(scenarioName.c_str(), frequencies[i], amplitudes[i], 24000);
    }
}


