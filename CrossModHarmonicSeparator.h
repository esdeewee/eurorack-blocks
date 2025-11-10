// CrossModHarmonicSeparator.h

#pragma once

#if defined (CROSSMOD_TEST_ENV_HEADER)
   #include CROSSMOD_TEST_ENV_HEADER
#else
   #include "artifacts/CrossModHarmonicSeparatorUi.h"
   //#include "artifacts/CrossModHarmonicSeparatorData.h" Uncomment if you use some 'data' resources
   #include "erb/erb.h"
#endif

#include "PitchDetector.h"


struct CrossModHarmonicSeparator
{
   // The UI elements defined in CrossModHarmonicSeparator.erbui are in 'ui'
   CrossModHarmonicSeparatorUi ui;

   // The data resources defined in CrossModHarmonicSeparator.erbb are in 'data'
   // Uncomment if you use some 'data' resources
   //CrossModHarmonicSeparatorData data;

   void  init ();
   void  process ();

   // Put here your DSP objects
   PitchDetector pitch_detector { 48000.0f };
   bool          has_detected_pitch = false;
   float         detected_pitch_hz = 0.0f;
};
