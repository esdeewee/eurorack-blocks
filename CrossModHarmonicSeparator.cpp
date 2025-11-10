// CrossModHarmonicSeparator.cpp

#include "CrossModHarmonicSeparator.h"

#include <cstddef>



void  CrossModHarmonicSeparator::init ()
{
   // This function is called once, before the first 'process' is called.
   // you can setup your module here.
   pitch_detector.reset();
   has_detected_pitch = false;
   detected_pitch_hz = 0.0f;
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
   }
   else
   {
      detected_pitch_hz = pitch_detector.getLastPitchHz();
   }
}
