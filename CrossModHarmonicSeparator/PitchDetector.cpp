#include "PitchDetector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace
{
constexpr float kDefaultThreshold      = 0.10f;
constexpr float kSilenceRmsThreshold   = 1.0e-5f;
constexpr float kMinFrequency          = 40.0f;
constexpr float kMaxFrequency          = 6000.0f;

float parabolicInterpolateDifference(const std::vector<float>& diff, std::size_t tau, std::size_t maxTau)
{
    if (tau == 0 || tau >= maxTau) {
        return static_cast<float>(tau);
    }
    if (tau + 1 > maxTau) {
        return static_cast<float>(tau);
    }

    const float prev = diff[tau - 1];
    const float cur  = diff[tau];
    const float next = diff[tau + 1];
    const float denominator = prev - 2.0f * cur + next;
    if (std::fabs(denominator) < 1.0e-12f) {
        return static_cast<float>(tau);
    }

    const float delta = 0.5f * (prev - next) / denominator;
    return static_cast<float>(tau) + delta;
}
} // namespace

PitchDetector::PitchDetector(float sampleRate, std::size_t analysisSize, std::size_t hopSize)
    : sample_rate_(sampleRate)
    , analysis_size_(analysisSize)
    , hop_size_(hopSize)
    , accumulator_()
    , difference_(analysis_size_, 0.0f)
    , cmnd_(analysis_size_, 0.0f)
{
    accumulator_.reserve(analysis_size_ * 2);
}

void PitchDetector::reset()
{
    accumulator_.clear();
    std::fill(difference_.begin(), difference_.end(), 0.0f);
    std::fill(cmnd_.begin(), cmnd_.end(), 0.0f);
    current_pitch_hz_.reset();
    last_valid_pitch_hz_.reset();
}

void PitchDetector::processBuffer(const float* buffer, std::size_t length)
{
    accumulator_.insert(accumulator_.end(), buffer, buffer + length);

    const std::size_t windowLength = std::min<std::size_t>(accumulator_.size(), analysis_size_);
    const std::size_t maxTauNeeded = static_cast<std::size_t>(sample_rate_ / kMinFrequency);
    const std::size_t minWindow = std::max<std::size_t>(maxTauNeeded, hop_size_ * 2);

    if (windowLength < minWindow) {
        current_pitch_hz_.reset();
        return;
    }

    const float* windowStart = accumulator_.data() + (accumulator_.size() - windowLength);
    const float rms = computeRms(windowStart, windowLength);

    if (rms < kSilenceRmsThreshold) {
        current_pitch_hz_.reset();
        if (accumulator_.size() > analysis_size_ + hop_size_) {
            accumulator_.erase(accumulator_.begin(), accumulator_.begin() + hop_size_);
        }
        return;
    }

    float pitchHz = 0.0f;
    if (detectPitch(windowStart, windowLength, pitchHz)) {
        current_pitch_hz_ = pitchHz;
        last_valid_pitch_hz_ = pitchHz;
    } else {
        current_pitch_hz_.reset();
    }

    if (accumulator_.size() > analysis_size_ + hop_size_) {
        accumulator_.erase(accumulator_.begin(), accumulator_.begin() + hop_size_);
    }
}

bool PitchDetector::hasPitch() const
{
    return current_pitch_hz_.has_value();
}

float PitchDetector::getCurrentPitchHz() const
{
    return current_pitch_hz_.value_or(0.0f);
}

float PitchDetector::getLastPitchHz() const
{
    return last_valid_pitch_hz_.value_or(0.0f);
}

bool PitchDetector::detectPitch(const float* windowStart, std::size_t windowLength, float& outPitchHz)
{
    if (windowLength < 3) {
        return false;
    }

    const std::size_t maxTauLimit = (windowLength > 1) ? (windowLength - 1) : 0;
    const std::size_t maxTau = std::min<std::size_t>(maxTauLimit, static_cast<std::size_t>(sample_rate_ / kMinFrequency));
    const std::size_t minTau = std::max<std::size_t>(2, static_cast<std::size_t>(sample_rate_ / kMaxFrequency));
    if (maxTau <= minTau) {
        return false;
    }

    std::fill(difference_.begin(), difference_.end(), 0.0f);
    std::fill(cmnd_.begin(), cmnd_.end(), 0.0f);

    for (std::size_t tau = 1; tau <= maxTau; ++tau) {
        float sum = 0.0f;
        const std::size_t limit = windowLength - tau;
        for (std::size_t i = 0; i < limit; ++i) {
            const float delta = windowStart[i] - windowStart[i + tau];
            sum += delta * delta;
        }
        difference_[tau] = sum;
    }

    float runningSum = 0.0f;
    cmnd_[0] = 1.0f;
    for (std::size_t tau = 1; tau <= maxTau; ++tau) {
        runningSum += difference_[tau];
        if (runningSum == 0.0f) {
            cmnd_[tau] = 1.0f;
        } else {
            cmnd_[tau] = difference_[tau] * static_cast<float>(tau) / runningSum;
        }
    }

    std::size_t bestTau = 0;
    bool candidateFound = false;
    for (std::size_t tau = minTau; tau <= maxTau; ++tau) {
        if (cmnd_[tau] < kDefaultThreshold) {
            while (tau + 1 <= maxTau && cmnd_[tau + 1] < cmnd_[tau]) {
                ++tau;
            }
            bestTau = tau;
            candidateFound = true;
            break;
        }
    }

    if (!candidateFound || bestTau == 0) {
        return false;
    }

    std::size_t refinedTauIndex = bestTau;
    while (refinedTauIndex / 2 >= minTau) {
        const std::size_t half = refinedTauIndex / 2;
        if (cmnd_[half] + 1.0e-3f < cmnd_[refinedTauIndex]) {
            refinedTauIndex = half;
        } else {
            break;
        }
    }

    const float refinedTau = parabolicInterpolateDifference(difference_, refinedTauIndex, maxTau);

    if (refinedTau <= 0.0f) {
        return false;
    }

    const float rawPitchHz = sample_rate_ / refinedTau;
    const float smoothingAlpha = 0.7f;
    if (last_valid_pitch_hz_) {
        outPitchHz = smoothingAlpha * rawPitchHz + (1.0f - smoothingAlpha) * last_valid_pitch_hz_.value();
    } else {
        outPitchHz = rawPitchHz;
    }
    if (!std::isfinite(outPitchHz) || outPitchHz <= 0.0f) {
        return false;
    }

    if (outPitchHz < kMinFrequency || outPitchHz > kMaxFrequency) {
        return false;
    }

    return true;
}

float PitchDetector::computeRms(const float* windowStart, std::size_t length) const
{
    double sumSquares = 0.0;
    for (std::size_t i = 0; i < length; ++i) {
        sumSquares += static_cast<double>(windowStart[i]) * static_cast<double>(windowStart[i]);
    }

    if (length == 0) {
        return 0.0f;
    }

    return static_cast<float>(std::sqrt(sumSquares / static_cast<double>(length)));
}


