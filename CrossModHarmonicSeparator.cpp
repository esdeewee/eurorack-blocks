// CrossModHarmonicSeparator.cpp

#include "CrossModHarmonicSeparator.h"

#include <array>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <limits>
#include <vector>
#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#endif

namespace
{
constexpr float kMinNormalizationValue = 1.0e-6f;

int64_t floorDiv(int64_t numerator, int64_t denominator)
{
   if (denominator <= 0)
   {
      return 0;
   }
   if (numerator >= 0)
   {
      return numerator / denominator;
   }
   return -(( -numerator + denominator - 1) / denominator);
}

int64_t ceilDiv(int64_t numerator, int64_t denominator)
{
   if (denominator <= 0)
   {
      return 0;
   }
   if (numerator >= 0)
   {
      return (numerator + denominator - 1) / denominator;
   }
   return -(( -numerator) / denominator);
}

std::vector<float> computeNormalizationPattern(const std::vector<float>& window,
                                               std::size_t hop);

#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
std::vector<float> computeWindowSumPattern(const std::vector<float>& window,
                                           std::size_t hop);
std::vector<int> computeOverlapCountPattern(const std::vector<float>& window,
                                            std::size_t hop);
#endif

std::vector<float> computeNormalizationPattern(const std::vector<float>& window,
                                               std::size_t hop)
{
   const std::size_t analysis_size = window.size();
   if (analysis_size == 0 || hop == 0)
   {
      return std::vector<float>(analysis_size == 0 ? 1 : analysis_size, 1.0f);
   }

   // Compute normalization pattern for overlap-add reconstruction
   // For each position n in the overlap buffer, calculate the sum of squared window values
   // from all overlapping frames that contribute to that position
   std::vector<float> pattern(analysis_size, 0.0f);
   
   const int64_t hop64 = static_cast<int64_t>(hop);
   const int64_t analysis64 = static_cast<int64_t>(analysis_size);
   
   // For each position in the overlap buffer
   for (std::size_t n = 0; n < analysis_size; ++n)
   {
      const int64_t n64 = static_cast<int64_t>(n);
      
      // Find all frames k that contribute to position n
      // Frame k starts at offset k * hop and covers positions [k*hop, k*hop + analysis_size)
      // Position n is covered by frame k if: k*hop <= n < k*hop + analysis_size
      // Rearranging: n - analysis_size < k*hop <= n
      // So: ceil((n - analysis_size + 1) / hop) <= k <= floor(n / hop)
      const int64_t min_k = ceilDiv(n64 - analysis64 + 1, hop64);
      const int64_t max_k = floorDiv(n64, hop64);

      float sum_sq = 0.0f;
      int count = 0;
      for (int64_t k = min_k; k <= max_k; ++k)
      {
         const int64_t frame_start = k * hop64;
         const int64_t pos_in_frame = n64 - frame_start;
         if (pos_in_frame >= 0 && pos_in_frame < analysis64)
         {
            const float w = window[static_cast<std::size_t>(pos_in_frame)];
            sum_sq += w * w;
            ++count;
         }
      }
      
      // In steady state, every position should have at least one contributing frame
      // If sum is too small, use window value as fallback (shouldn't happen in steady state)
      if (sum_sq < kMinNormalizationValue)
      {
         if (n < window.size())
         {
            const float w = window[n];
            sum_sq = std::max(w * w, kMinNormalizationValue);
         }
         else
         {
            sum_sq = kMinNormalizationValue;
         }
      }
      pattern[n] = sum_sq;
   }

   return pattern;
}

      #if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
// Compute sum of windows (not squared) for COLA verification
std::vector<float> computeWindowSumPattern(const std::vector<float>& window,
                                           std::size_t hop)
{
   const std::size_t analysis_size = window.size();
   if (analysis_size == 0 || hop == 0)
   {
      return std::vector<float>(analysis_size == 0 ? 1 : analysis_size, 1.0f);
   }

   std::vector<float> pattern(analysis_size, 0.0f);
   const int64_t hop64 = static_cast<int64_t>(hop);
   const int64_t analysis64 = static_cast<int64_t>(analysis_size);
   
   for (std::size_t n = 0; n < analysis_size; ++n)
   {
      const int64_t n64 = static_cast<int64_t>(n);
      const int64_t min_k = ceilDiv(n64 - analysis64 + 1, hop64);
      const int64_t max_k = floorDiv(n64, hop64);

      float sum = 0.0f;
      for (int64_t k = min_k; k <= max_k; ++k)
      {
         const int64_t frame_start = k * hop64;
         const int64_t pos_in_frame = n64 - frame_start;
         if (pos_in_frame >= 0 && pos_in_frame < analysis64)
         {
            const float w = window[static_cast<std::size_t>(pos_in_frame)];
            sum += w;
         }
      }
      pattern[n] = sum;
   }

   return pattern;
}

// Compute number of overlapping frames at each position
std::vector<int> computeOverlapCountPattern(const std::vector<float>& window,
                                            std::size_t hop)
{
   const std::size_t analysis_size = window.size();
   if (analysis_size == 0 || hop == 0)
   {
      return std::vector<int>(analysis_size == 0 ? 1 : analysis_size, 1);
   }

   std::vector<int> pattern(analysis_size, 0);
   const int64_t hop64 = static_cast<int64_t>(hop);
   const int64_t analysis64 = static_cast<int64_t>(analysis_size);
   
   for (std::size_t n = 0; n < analysis_size; ++n)
   {
      const int64_t n64 = static_cast<int64_t>(n);
      const int64_t min_k = ceilDiv(n64 - analysis64 + 1, hop64);
      const int64_t max_k = floorDiv(n64, hop64);

      int count = 0;
      for (int64_t k = min_k; k <= max_k; ++k)
      {
         const int64_t frame_start = k * hop64;
         const int64_t pos_in_frame = n64 - frame_start;
         if (pos_in_frame >= 0 && pos_in_frame < analysis64)
         {
            ++count;
         }
      }
      pattern[n] = count;
   }

   return pattern;
}
#endif
} // namespace



void  CrossModHarmonicSeparator::init ()
{
   // This function is called once, before the first 'process' is called.
   // you can setup your module here.
   pitch_detector.reset();
   harmonic_series.reset();
   spectral_separator.reset();
   has_detected_pitch = false;
   detected_pitch_hz = 0.0f;
   has_harmonic_series = false;
   fundamental_group.clear();
   prime_group.clear();
   composite_group.clear();
   prime_numbers.clear();
   composite_numbers.clear();
   has_spectral_separation = false;
   fundamental_time.clear();
   prime_time.clear();
   composite_time.clear();
   energy_total = 0.0f;
   energy_fundamental = 0.0f;
   energy_prime = 0.0f;
   energy_composite = 0.0f;
   reconstruction_error = 0.0f;
   const std::size_t analysis_size = spectral_separator.analysisSize();
   fundamental_overlap.assign(analysis_size, 0.0f);
   prime_overlap.assign(analysis_size, 0.0f);
   composite_overlap.assign(analysis_size, 0.0f);
   normalization_overlap.assign(analysis_size, 0.0f);
   normalization_pattern.assign(analysis_size, 0.0f);

   const auto& window = spectral_separator.window();
   const std::size_t hop = spectral_separator.hopSize();
   if (analysis_size > 0)
   {
      if (hop > 0 && analysis_size > hop)
      {
         normalization_pattern = computeNormalizationPattern(window, hop);
         #if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
         // Also compute sum of windows (not squared) for COLA verification
         debug_window_sum_pattern = computeWindowSumPattern(window, hop);
         debug_overlap_count_pattern = computeOverlapCountPattern(window, hop);
         #endif
      }
      else
      {
         for (std::size_t n = 0; n < analysis_size; ++n)
         {
            const float w = window[n];
            normalization_pattern[n] = std::max(w * w, kMinNormalizationValue);
         }
         #if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
         debug_window_sum_pattern.assign(analysis_size, 0.0f);
         debug_overlap_count_pattern.assign(analysis_size, 0);
         for (std::size_t n = 0; n < analysis_size; ++n)
         {
            debug_window_sum_pattern[n] = window[n];
            debug_overlap_count_pattern[n] = 1;
         }
         #endif
      }
   }
   else
   {
      normalization_pattern.assign(1, 1.0f);
      #if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
      debug_window_sum_pattern.assign(1, 1.0f);
      debug_overlap_count_pattern.assign(1, 1);
      #endif
   }
   normalization_overlap.assign(analysis_size, 0.0f);

   fundamental_processor.init(48000.0f);
   fundamental_processor.setLevel(fundamental_params.level);
   fundamental_processor.setTone(fundamental_params.tone);
   fundamental_processor.setSpread(fundamental_params.spread);
   fundamental_processor.reset();
   fundamental_frame.fill(0.0f);
   fundamental_processed_block.fill(0.0f);

   prime_processor.init(48000.0f);
   prime_processor.setExciter(prime_params.exciter);
   prime_processor.setBrightness(prime_params.brightness);
   prime_processor.setPhaseRotation(prime_params.phase_rotation);
   prime_processor.reset();
   prime_frame.fill(0.0f);
   prime_processed_block.fill(0.0f);

   composite_processor.init(48000.0f);
   composite_processor.setSaturation(composite_params.saturation);
   composite_processor.setWarmth(composite_params.warmth);
   composite_processor.setDensity(composite_params.density);
   composite_processor.reset();
   composite_frame.fill(0.0f);
   composite_processed_block.fill(0.0f);
}


void  CrossModHarmonicSeparator::process ()
{
   // This function is called regularly every buffer size
   // get your audio input(s) if any, and write to your audio output(s)

   std::array<float, erb_BUFFER_SIZE> input_block {};
   fundamental_frame.fill(0.0f);

   for (std::size_t i = 0 ; i < erb_BUFFER_SIZE ; ++i)
   {
      const float sample = ui.audio_in [i];
      ui.audio_out [i] = sample;
      input_block [i] = sample;
   }

   pitch_detector.processBuffer(input_block.data(), erb_BUFFER_SIZE);
   has_detected_pitch = pitch_detector.hasPitch();

   const std::size_t analysis_size = spectral_separator.analysisSize();
   const std::size_t hop = std::min<std::size_t>(spectral_separator.hopSize(), erb_BUFFER_SIZE);

   if (has_detected_pitch)
   {
      detected_pitch_hz = pitch_detector.getCurrentPitchHz();
      harmonic_series.update(detected_pitch_hz);
      has_harmonic_series = harmonic_series.hasSeries();
   }
   else
   {
      detected_pitch_hz = pitch_detector.getLastPitchHz();
      harmonic_series.reset();
      has_harmonic_series = false;
   }

   if (has_harmonic_series)
   {
      const auto& fundamental = harmonic_series.fundamentalFrequencies();
      const auto& primes = harmonic_series.primeFrequencies();
      const auto& composites = harmonic_series.compositeFrequencies();
      const auto& primeNums = harmonic_series.primeNumbers();
      const auto& compositeNums = harmonic_series.compositeNumbers();

      fundamental_group.assign(fundamental.begin(), fundamental.end());
      prime_group.assign(primes.begin(), primes.end());
      composite_group.assign(composites.begin(), composites.end());
      prime_numbers.assign(primeNums.begin(), primeNums.end());
      composite_numbers.assign(compositeNums.begin(), compositeNums.end());

      spectral_separator.setHarmonicData(detected_pitch_hz, prime_numbers, composite_numbers);
      if (spectral_separator.processBuffer(input_block.data(), erb_BUFFER_SIZE) && spectral_separator.hasResult())
      {
         has_spectral_separation = true;
         const auto& fundamental_domain = spectral_separator.fundamentalTimeDomain();
         const auto& prime_domain = spectral_separator.primeTimeDomain();
         const auto& composite_domain = spectral_separator.compositeTimeDomain();

         fundamental_time.assign(fundamental_domain.begin(), fundamental_domain.end());
         prime_time.assign(prime_domain.begin(), prime_domain.end());
         composite_time.assign(composite_domain.begin(), composite_domain.end());

         const auto& window = spectral_separator.window();
         if (fundamental_overlap.size() != analysis_size)
         {
            fundamental_overlap.assign(analysis_size, 0.0f);
            prime_overlap.assign(analysis_size, 0.0f);
            composite_overlap.assign(analysis_size, 0.0f);
            normalization_overlap.assign(analysis_size, 0.0f);
            normalization_pattern.assign(analysis_size, 0.0f);

            const std::size_t hop = spectral_separator.hopSize();
            if (hop > 0 && analysis_size > hop)
            {
               normalization_pattern = computeNormalizationPattern(window, hop);
               #if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
               // Also compute sum of windows (not squared) for COLA verification
               debug_window_sum_pattern = computeWindowSumPattern(window, hop);
               debug_overlap_count_pattern = computeOverlapCountPattern(window, hop);
               #endif
            }
            else
            {
               for (std::size_t n = 0; n < analysis_size; ++n)
               {
                  const float w = window[n];
                  normalization_pattern[n] = std::max(w * w, kMinNormalizationValue);
               }
               #if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
               debug_window_sum_pattern.assign(analysis_size, 0.0f);
               debug_overlap_count_pattern.assign(analysis_size, 0);
               for (std::size_t n = 0; n < analysis_size; ++n)
               {
                  debug_window_sum_pattern[n] = window[n];
                  debug_overlap_count_pattern[n] = 1;
               }
               #endif
            }
            normalization_overlap.assign(analysis_size, 0.0f);
         }

        energy_total = spectral_separator.totalEnergy();
        energy_fundamental = spectral_separator.fundamentalEnergy();
        energy_prime = spectral_separator.primeEnergy();
        energy_composite = spectral_separator.compositeEnergy();
        reconstruction_error = spectral_separator.reconstructionError();
#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
        debug_total_separated_energy = static_cast<double>(energy_fundamental + energy_prime + energy_composite);
        debug_energy_total_latest = energy_total;
        debug_spectral_ratio = (energy_total > 0.0f && debug_total_separated_energy > 0.0)
           ? (debug_total_separated_energy / energy_total)
           : 0.0;
#endif
      }
      else
      {
         has_spectral_separation = false;
      }
   }
   else
   {
      fundamental_group.clear();
      prime_group.clear();
      composite_group.clear();
      prime_numbers.clear();
      composite_numbers.clear();
      spectral_separator.reset();
      has_spectral_separation = false;
      fundamental_time.clear();
      prime_time.clear();
      composite_time.clear();
        energy_total = 0.0f;
        energy_fundamental = 0.0f;
        energy_prime = 0.0f;
        energy_composite = 0.0f;
        reconstruction_error = 0.0f;
#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
        debug_total_separated_energy = 0.0;
        debug_energy_total_latest = 0.0;
        debug_spectral_ratio = 0.0;
#endif
      std::fill(fundamental_overlap.begin(), fundamental_overlap.end(), 0.0f);
      std::fill(prime_overlap.begin(), prime_overlap.end(), 0.0f);
      std::fill(composite_overlap.begin(), composite_overlap.end(), 0.0f);
      normalization_overlap.assign(analysis_size, 0.0f);
   }

   
   double hop_input_energy_measure = 0.0;
#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
   // Measure energy before normalization
   debug_pre_normalization_energy = 0.0;
   debug_normalization_values.clear();
   debug_normalization_values.reserve(hop);
#endif
   for (std::size_t i = 0; i < hop; ++i)
   {
#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
      if (i < fundamental_overlap.size())
      {
         const double val = static_cast<double>(fundamental_overlap[i]);
         debug_pre_normalization_energy += val * val;
      }
      if (i < prime_overlap.size())
      {
         const double val = static_cast<double>(prime_overlap[i]);
         debug_pre_normalization_energy += val * val;
      }
      if (i < composite_overlap.size())
      {
         const double val = static_cast<double>(composite_overlap[i]);
         debug_pre_normalization_energy += val * val;
      }
#endif
      if (i < erb_BUFFER_SIZE)
      {
         const double val = static_cast<double>(ui.audio_in[i]);
         hop_input_energy_measure += val * val;
      }
   }
#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
   debug_input_energy = hop_input_energy_measure;
#endif
   
   // Apply normalization for overlap-add reconstruction
   // The normalization pattern contains the sum of squared window values at each position
   // We divide by this to compensate for the energy scaling from overlapping windows
#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
   double norm_min_measure = std::numeric_limits<double>::max();
   double norm_max_measure = 0.0;
   double norm_sum_measure = 0.0;
#endif
   for (std::size_t i = 0 ; i < hop ; ++i)
   {
      float norm = (i < normalization_overlap.size()) ? normalization_overlap[i] : 0.0f;
      if (norm <= kMinNormalizationValue)
      {
         norm = 1.0f;
      }
      #if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
      debug_normalization_values.push_back(static_cast<double>(norm));
      norm_min_measure = std::min(norm_min_measure, static_cast<double>(norm));
      norm_max_measure = std::max(norm_max_measure, static_cast<double>(norm));
      norm_sum_measure += norm;
      #endif
      
      const float inv_norm = 1.0f / norm;
      float fund_sample = (i < fundamental_overlap.size()) ? fundamental_overlap[i] * inv_norm : 0.0f;
      float prime_sample = (i < prime_overlap.size()) ? prime_overlap[i] * inv_norm : 0.0f;
      float comp_sample = (i < composite_overlap.size()) ? composite_overlap[i] * inv_norm : 0.0f;

      fund_sample *= CrossModHarmonicSeparator::kEnergyNormalizationGain;
      prime_sample *= CrossModHarmonicSeparator::kEnergyNormalizationGain;
      comp_sample *= CrossModHarmonicSeparator::kEnergyNormalizationGain;

      // Perform normalized overlap-add using the live sum of squared windows stored in normalization_overlap.
      ui.fundamental_debug [i] = fund_sample;
      if (i < fundamental_frame.size())
      {
         fundamental_frame[i] = fund_sample;
      }
      ui.prime_debug [i] = prime_sample;
      if (i < prime_frame.size())
      {
         prime_frame[i] = prime_sample;
      }
      ui.composite_debug [i] = comp_sample;
      if (i < composite_frame.size())
      {
         composite_frame[i] = comp_sample;
      }
   }
   
   double hop_output_energy_measure = 0.0;
   for (std::size_t i = 0; i < hop; ++i)
   {
      const double val_fund = static_cast<double>(ui.fundamental_debug[i]);
      const double val_prime = static_cast<double>(ui.prime_debug[i]);
      const double val_comp = static_cast<double>(ui.composite_debug[i]);
      hop_output_energy_measure += val_fund * val_fund + val_prime * val_prime + val_comp * val_comp;
   }
#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
   debug_post_normalization_energy = hop_output_energy_measure;
#endif
   for (std::size_t i = hop ; i < erb_BUFFER_SIZE ; ++i)
   {
      ui.fundamental_debug [i] = 0.0f;
      ui.prime_debug [i] = 0.0f;
      ui.composite_debug [i] = 0.0f;
      fundamental_frame[i] = 0.0f;
      prime_frame[i] = 0.0f;
      composite_frame[i] = 0.0f;
   }

   fundamental_processor.processBuffer(fundamental_frame.data(),
                                       fundamental_processed_block.data(),
                                       fundamental_frame.size());

   prime_processor.processBuffer(prime_frame.data(),
                                 prime_processed_block.data(),
                                 prime_frame.size());
   composite_processor.processBuffer(composite_frame.data(),
                                     composite_processed_block.data(),
                                     composite_frame.size());

#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
   if (hop > 0)
   {
      debug_norm_min = (norm_min_measure == std::numeric_limits<double>::max()) ? 0.0 : norm_min_measure;
      debug_norm_max = norm_max_measure;
      debug_norm_avg = norm_sum_measure / static_cast<double>(hop);
   }
   debug_last_normalization_overlap = normalization_overlap;
   debug_last_fundamental_overlap = fundamental_overlap;
   debug_last_prime_overlap = prime_overlap;
   debug_last_composite_overlap = composite_overlap;
   debug_last_input_buffer.assign(ui.audio_in.begin(), ui.audio_in.end());
#endif

   auto shiftBuffer = [hop](std::vector<float>& buffer) {
      if (buffer.empty())
      {
         return;
      }
      if (hop >= buffer.size())
      {
         std::fill(buffer.begin(), buffer.end(), 0.0f);
         return;
      }
      // Shift left by hop samples (source is to the right of destination, so std::copy is safe)
      std::copy(buffer.begin() + hop, buffer.end(), buffer.begin());
      std::fill(buffer.end() - hop, buffer.end(), 0.0f);
   };

   shiftBuffer(fundamental_overlap);
   shiftBuffer(prime_overlap);
   shiftBuffer(composite_overlap);
   shiftBuffer(normalization_overlap);

   // Add new windowed frames to overlap buffers after shifting
   // After shifting, position 0 is where the new frame should start (overlapping with remaining tail)
   if (has_spectral_separation)
   {
      const auto& fundamental_domain = spectral_separator.fundamentalTimeDomain();
      const auto& prime_domain = spectral_separator.primeTimeDomain();
      const auto& composite_domain = spectral_separator.compositeTimeDomain();
      const auto& window = spectral_separator.window();
      
      // Add windowed frame starting at position 0 (overlaps with tail from previous frames)
      // Note: After shift, positions 0..analysis_size-hop-1 contain tail from previous frames
      // Positions analysis_size-hop..analysis_size-1 are zero (were shifted out)
      // We add new frame starting at position 0, which overlaps with the tail
      //
      // FUNDAMENTAL FIX: The IFFT output has 1/N normalization applied in computeInverse.
      // This causes energy to be scaled by 1/N^2. The cleanest solution is to adjust
      // the normalization pattern to account for this, rather than multiplying here.
      // We'll handle the compensation in the normalization step instead.
      for (std::size_t n = 0; n < analysis_size; ++n)
      {
         const float win = window[n];
         // Keep IFFT output as-is (with 1/N normalization)
         // The normalization pattern will be adjusted to account for 1/N^2 scaling
         const float scaled = fundamental_domain[n] * win;
         const float scaled_prime = prime_domain[n] * win;
         const float scaled_composite = composite_domain[n] * win;
         if (n < fundamental_overlap.size())
         {
            fundamental_overlap[n] += scaled;
         }
         if (n < prime_overlap.size())
         {
            prime_overlap[n] += scaled_prime;
         }
         if (n < composite_overlap.size())
         {
            composite_overlap[n] += scaled_composite;
         }
         if (n < normalization_overlap.size())
         {
            normalization_overlap[n] += win * win;
         }
      }
   }
}

#if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
void CrossModHarmonicSeparator::dumpDebugState(const char* label) const
{
   std::ofstream log_stream("crossmod_debug_state.log", std::ios::app);
   auto emit = [&](const std::string& text)
   {
      std::cout << text << '\n';
      if (log_stream.is_open())
      {
         log_stream << text << '\n';
      }
   };
   auto formatFloatSlice = [](const char* name, const std::vector<float>& data, std::size_t limit = 16)
   {
      std::ostringstream oss;
      oss << name << " (size=" << data.size()
          << ", first " << std::min(limit, data.size()) << "): ";
      const std::size_t count = std::min(limit, data.size());
      for (std::size_t i = 0; i < count; ++i)
      {
         if (i > 0) oss << ", ";
         oss << data[i];
      }
      if (data.size() > count)
      {
         oss << ", ...";
      }
      return oss.str();
   };
   auto formatDoubleSlice = [](const char* name, const std::vector<double>& data, std::size_t limit = 16)
   {
      std::ostringstream oss;
      oss << name << " (size=" << data.size()
          << ", first " << std::min(limit, data.size()) << "): ";
      const std::size_t count = std::min(limit, data.size());
      for (std::size_t i = 0; i < count; ++i)
      {
         if (i > 0) oss << ", ";
         oss << data[i];
      }
      if (data.size() > count)
      {
         oss << ", ...";
      }
      return oss.str();
   };

   std::ostringstream header;
   header << "\n===== CrossMod Debug State [" << (label ? label : "unnamed") << "] =====";
   emit(header.str());
   emit((std::ostringstream() << std::fixed << std::setprecision(6)
        << "Input energy (hop): " << debug_input_energy).str());
   emit((std::ostringstream() << "Pre-normalization energy: " << debug_pre_normalization_energy).str());
   emit((std::ostringstream() << "Post-normalization energy: " << debug_post_normalization_energy).str());
   emit((std::ostringstream() << "Spectral ratio: " << debug_spectral_ratio
         << " (total separated=" << debug_total_separated_energy
         << ", energy_total=" << debug_energy_total_latest << ")").str());
   emit((std::ostringstream() << "Normalization stats: min=" << debug_norm_min
         << ", max=" << debug_norm_max
         << ", avg=" << debug_norm_avg).str());

   emit(formatFloatSlice("Normalization overlap snapshot", debug_last_normalization_overlap));
   emit(formatFloatSlice("Fundamental overlap snapshot", debug_last_fundamental_overlap));
   emit(formatFloatSlice("Prime overlap snapshot", debug_last_prime_overlap));
   emit(formatFloatSlice("Composite overlap snapshot", debug_last_composite_overlap));
   emit(formatFloatSlice("Input buffer snapshot", debug_last_input_buffer));
   emit(formatDoubleSlice("Recent normalization values", debug_normalization_values));
   emit("=============================================");
}
#endif

void CrossModHarmonicSeparator::setFundamentalLevel(float value)
{
   const float clamped = std::clamp(value, 0.0f, 1.0f);
   fundamental_params.level = clamped;
   fundamental_processor.setLevel(clamped);
}

void CrossModHarmonicSeparator::setFundamentalTone(float value)
{
   const float clamped = std::clamp(value, 0.0f, 1.0f);
   fundamental_params.tone = clamped;
   fundamental_processor.setTone(clamped);
}

void CrossModHarmonicSeparator::setFundamentalSpread(float value)
{
   const float clamped = std::clamp(value, 0.0f, 1.0f);
   fundamental_params.spread = clamped;
   fundamental_processor.setSpread(clamped);
}

void CrossModHarmonicSeparator::setPrimeExciter(float value)
{
   const float clamped = std::clamp(value, 0.0f, 1.0f);
   prime_params.exciter = clamped;
   prime_processor.setExciter(clamped);
}

void CrossModHarmonicSeparator::setPrimeBrightness(float value)
{
   const float clamped = std::clamp(value, 0.0f, 1.0f);
   prime_params.brightness = clamped;
   prime_processor.setBrightness(clamped);
}

void CrossModHarmonicSeparator::setPrimePhaseRotation(float value)
{
   const float clamped = std::clamp(value, 0.0f, 1.0f);
   prime_params.phase_rotation = clamped;
   prime_processor.setPhaseRotation(clamped);
}

void CrossModHarmonicSeparator::setCompositeSaturation(float value)
{
   const float clamped = std::clamp(value, 0.0f, 1.0f);
   composite_params.saturation = clamped;
   composite_processor.setSaturation(clamped);
}

void CrossModHarmonicSeparator::setCompositeWarmth(float value)
{
   const float clamped = std::clamp(value, 0.0f, 1.0f);
   composite_params.warmth = clamped;
   composite_processor.setWarmth(clamped);
}

void CrossModHarmonicSeparator::setCompositeDensity(float value)
{
   const float clamped = std::clamp(value, 0.0f, 1.0f);
   composite_params.density = clamped;
   composite_processor.setDensity(clamped);
}