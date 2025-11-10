#pragma once

#include <array>
#include <cstddef>

constexpr std::size_t erb_BUFFER_SIZE = 48;
constexpr std::size_t erb_SAMPLE_RATE = 48000;

struct CrossModHarmonicSeparatorUi
{
    std::array<float, erb_BUFFER_SIZE> audio_in {};
    std::array<float, erb_BUFFER_SIZE> audio_out {};
};

namespace erb
{
using Buffer = std::array<float, erb_BUFFER_SIZE>;
}


