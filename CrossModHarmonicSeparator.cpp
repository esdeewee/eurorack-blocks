// CrossModHarmonicSeparator.cpp

#include "CrossModHarmonicSeparator.h"

#include <cstddef>
#include <array>
#include <algorithm>



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
}


void  CrossModHarmonicSeparator::process ()
{
   // This function is called regularly every buffer size
   // get your audio input(s) if any, and write to your audio output(s)

   std::array<float, erb_BUFFER_SIZE> input_block {};

   for (std::size_t i = 0 ; i < erb_BUFFER_SIZE ; ++i)
   {
      const float sample = ui.audio_in [i];
      ui.audio_out [i] = sample;
      input_block [i] = sample;
   }

   pitch_detector.processBuffer(input_block.data(), erb_BUFFER_SIZE);
   has_detected_pitch = pitch_detector.hasPitch();

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
         const std::size_t analysis_size = spectral_separator.analysisSize();
         if (fundamental_overlap.size() != analysis_size)
         {
            fundamental_overlap.assign(analysis_size, 0.0f);
            prime_overlap.assign(analysis_size, 0.0f);
            composite_overlap.assign(analysis_size, 0.0f);
            normalization_overlap.assign(analysis_size, 0.0f);
         }

         for (std::size_t n = 0; n < analysis_size; ++n)
         {
            const float win = window[n];
            fundamental_overlap[n] += fundamental_domain[n] * win;
            prime_overlap[n] += prime_domain[n] * win;
            composite_overlap[n] += composite_domain[n] * win;
            normalization_overlap[n] += win * win;
         }

         energy_total = spectral_separator.totalEnergy();
         energy_fundamental = spectral_separator.fundamentalEnergy();
         energy_prime = spectral_separator.primeEnergy();
         energy_composite = spectral_separator.compositeEnergy();
         reconstruction_error = spectral_separator.reconstructionError();
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
      std::fill(fundamental_overlap.begin(), fundamental_overlap.end(), 0.0f);
      std::fill(prime_overlap.begin(), prime_overlap.end(), 0.0f);
      std::fill(composite_overlap.begin(), composite_overlap.end(), 0.0f);
      std::fill(normalization_overlap.begin(), normalization_overlap.end(), 0.0f);
   }

   const std::size_t hop = std::min<std::size_t>(spectral_separator.hopSize(), erb_BUFFER_SIZE);
   for (std::size_t i = 0 ; i < hop ; ++i)
   {
      const float norm = (i < normalization_overlap.size() && normalization_overlap[i] > 1e-6f) ? normalization_overlap[i] : 1.0f;
      ui.fundamental_debug [i] = (i < fundamental_overlap.size()) ? fundamental_overlap[i] / norm : 0.0f;
      ui.prime_debug [i] = (i < prime_overlap.size()) ? prime_overlap[i] / norm : 0.0f;
      ui.composite_debug [i] = (i < composite_overlap.size()) ? composite_overlap[i] / norm : 0.0f;
   }
   for (std::size_t i = hop ; i < erb_BUFFER_SIZE ; ++i)
   {
      ui.fundamental_debug [i] = 0.0f;
      ui.prime_debug [i] = 0.0f;
      ui.composite_debug [i] = 0.0f;
   }

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
      std::move(buffer.begin() + hop, buffer.end(), buffer.begin());
      std::fill(buffer.end() - hop, buffer.end(), 0.0f);
   };

   shiftBuffer(fundamental_overlap);
   shiftBuffer(prime_overlap);
   shiftBuffer(composite_overlap);
   shiftBuffer(normalization_overlap);
}
