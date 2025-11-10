#include "PitchDetector.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

namespace
{
constexpr float kDefaultThreshold = 0.10f;
constexpr float kSilenceRmsThreshold = 1.0e-4f;
constexpr float kMinFrequency = 50.0f;
constexpr float kMaxFrequency = 2000.0f;

float parabolicInterpolate(float tau, float prev, float cur, float next)
{
    const float denominator = 2.0f * (prev - 2.0f * cur + next);
    if (denominator == 0.0f) {
        return tau;
    }

    const float delta = (prev - next) / denominator;
    return tau + delta;
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

    if (accumulator_.size() < analysis_size_) {
        return;
    }

    const float* windowStart = accumulator_.data() + (accumulator_.size() - analysis_size_);
    const float rms = computeRms(windowStart, analysis_size_);

    if (rms < kSilenceRmsThreshold) {
        current_pitch_hz_.reset();
        if (accumulator_.size() > analysis_size_ + hop_size_) {
            accumulator_.erase(accumulator_.begin(), accumulator_.begin() + hop_size_);
        }
        return;
    }

    float pitchHz = 0.0f;
    if (detectPitch(windowStart, pitchHz)) {
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

bool PitchDetector::detectPitch(const float* windowStart, float& outPitchHz)
{
    const std::size_t maxTau = std::min<std::size_t>(analysis_size_ / 2, static_cast<std::size_t>(sample_rate_ / kMinFrequency));
    const std::size_t minTau = std::max<std::size_t>(2, static_cast<std::size_t>(sample_rate_ / kMaxFrequency));
    const std::size_t windowLimit = analysis_size_ - 1;

    std::fill(difference_.begin(), difference_.end(), 0.0f);
    std::fill(cmnd_.begin(), cmnd_.end(), 0.0f);

    for (std::size_t tau = 1; tau <= maxTau; ++tau) {
        float sum = 0.0f;
        for (std::size_t i = 0; i < windowLimit - tau; ++i) {
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

    const float prev = cmnd_[bestTau > 1 ? bestTau - 1 : bestTau];
    const float cur = cmnd_[bestTau];
    const float next = cmnd_[bestTau + 1 <= maxTau ? bestTau + 1 : bestTau];
    const float refinedTau = parabolicInterpolate(static_cast<float>(bestTau), prev, cur, next);

    if (refinedTau <= 0.0f) {
        return false;
    }

    outPitchHz = sample_rate_ / refinedTau;
    return std::isfinite(outPitchHz) && outPitchHz > 0.0f;
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


