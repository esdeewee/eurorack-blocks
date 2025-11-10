// CrossModHarmonicSeparator.cpp

#include "CrossModHarmonicSeparator.h"

#include <cstddef>



void  CrossModHarmonicSeparator::init ()
{
   // This function is called once, before the first 'process' is called.
   // you can setup your module here.
   pitch_detector.reset();
   harmonic_series.reset();
   has_detected_pitch = false;
   detected_pitch_hz = 0.0f;
   has_harmonic_series = false;
   fundamental_group.clear();
   prime_group.clear();
   composite_group.clear();
   prime_numbers.clear();
   composite_numbers.clear();
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
   }
   else
   {
      fundamental_group.clear();
      prime_group.clear();
      composite_group.clear();
      prime_numbers.clear();
      composite_numbers.clear();
   }
}
