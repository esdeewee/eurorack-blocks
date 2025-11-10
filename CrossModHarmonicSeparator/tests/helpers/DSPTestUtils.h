#pragma once

#include <vector>
#include <cmath>
#include <numeric>

inline double computeRms(const std::vector<float> &buffer)
{
    if (buffer.empty()) {
        return 0.0;
    }

    double sumSquares = 0.0;
    for (float sample : buffer) {
        sumSquares += static_cast<double>(sample) * static_cast<double>(sample);
    }

    return std::sqrt(sumSquares / static_cast<double>(buffer.size()));
}

inline double computeEnergy(const std::vector<float>& buffer)
{
    double sum = 0.0;
    for (float sample : buffer) {
        sum += static_cast<double>(sample) * static_cast<double>(sample);
    }
    return sum;
}

inline double computeRmsDifference(const std::vector<float>& a,
                                   const std::vector<float>& b)
{
    if (a.size() != b.size() || a.empty()) {
        return 0.0;
    }

    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        const double diff = static_cast<double>(a[i]) - static_cast<double>(b[i]);
        sum += diff * diff;
    }
    return std::sqrt(sum / static_cast<double>(a.size()));
}


