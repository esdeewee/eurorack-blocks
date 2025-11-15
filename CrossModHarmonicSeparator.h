// CrossModHarmonicSeparator.h

#pragma once

#include <cstddef>
#include <vector>

#if defined (CROSSMOD_TEST_ENV_HEADER)
   #include CROSSMOD_TEST_ENV_HEADER
#else
#include "artifacts/CrossModHarmonicSeparatorUi.h"
//#include "artifacts/CrossModHarmonicSeparatorData.h" Uncomment if you use some 'data' resources
#include "erb/erb.h"
#endif

#include "CrossModHarmonicSeparator/PitchDetector.h"
#include "CrossModHarmonicSeparator/HarmonicSeriesGenerator.h"
#include "CrossModHarmonicSeparator/HarmonicSpectralSeparator.h"


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
   PitchDetector             pitch_detector { 48000.0f };
   HarmonicSeriesGenerator   harmonic_series { 48000.0f };
   HarmonicSpectralSeparator spectral_separator { 48000.0f };

   bool   has_detected_pitch = false;
   float  detected_pitch_hz = 0.0f;

   bool   has_harmonic_series = false;
   std::vector<float> fundamental_group;
   std::vector<float> prime_group;
   std::vector<float> composite_group;
   std::vector<std::size_t> prime_numbers;
   std::vector<std::size_t> composite_numbers;

   bool   has_spectral_separation = false;
   std::vector<float> fundamental_time;
   std::vector<float> prime_time;
   std::vector<float> composite_time;
   float energy_total = 0.0f;
   float energy_fundamental = 0.0f;
   float energy_prime = 0.0f;
   float energy_composite = 0.0f;
   float reconstruction_error = 0.0f;
   std::vector<float> fundamental_overlap;
   std::vector<float> prime_overlap;
   std::vector<float> composite_overlap;
   std::vector<float> normalization_overlap;
   std::vector<float> normalization_pattern;
   std::size_t normalization_phase = 0;
   
   // Debug members for energy analysis (only in test environment)
   #if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
   std::vector<float> debug_window_sum_pattern;  // Sum of windows (not squared) for COLA verification
   std::vector<int> debug_overlap_count_pattern;  // Number of overlapping frames at each position
   double debug_pre_normalization_energy = 0.0;   // Energy before normalization
   double debug_post_normalization_energy = 0.0;   // Energy after normalization
   double debug_input_energy = 0.0;                // Input energy for current buffer
   std::vector<double> debug_normalization_values; // Normalization values used for output
   #endif
};
