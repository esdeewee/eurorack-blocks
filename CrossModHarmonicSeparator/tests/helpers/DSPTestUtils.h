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


