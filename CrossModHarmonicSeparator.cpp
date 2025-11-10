// CrossModHarmonicSeparator.cpp

#include "CrossModHarmonicSeparator.h"

#include <cstddef>



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
}


void  CrossModHarmonicSeparator::process ()
{
   // This function is called regularly every buffer size
   // get your audio input(s) if any, and write to your audio output(s)

   for (std::size_t i = 0 ; i < erb_BUFFER_SIZE ; ++i)
   {
      ui.audio_out [i] = ui.audio_in [i];
   }

   pitch_detector.processBuffer(ui.audio_in.data(), erb_BUFFER_SIZE);
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
      if (spectral_separator.processBuffer(ui.audio_in.data(), erb_BUFFER_SIZE) && spectral_separator.hasResult())
      {
         has_spectral_separation = true;
         fundamental_time.assign(spectral_separator.fundamentalTimeDomain().begin(),
                                 spectral_separator.fundamentalTimeDomain().end());
         prime_time.assign(spectral_separator.primeTimeDomain().begin(),
                           spectral_separator.primeTimeDomain().end());
         composite_time.assign(spectral_separator.compositeTimeDomain().begin(),
                               spectral_separator.compositeTimeDomain().end());
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
   }
}
