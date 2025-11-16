// CrossModHarmonicSeparator.h

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#if defined (CROSSMOD_TEST_ENV_HEADER)
   #include CROSSMOD_TEST_ENV_HEADER
#else
#include "artifacts/CrossModHarmonicSeparatorUi.h"
//#include "artifacts/CrossModHarmonicSeparatorData.h" Uncomment if you use some 'data' resources
#include "erb/erb.h"
#endif

#include "CrossModHarmonicSeparator/FundamentalProcessor.h"
#include "CrossModHarmonicSeparator/PrimeProcessor.h"
#include "CrossModHarmonicSeparator/CompositeProcessor.h"
#include "CrossModHarmonicSeparator/PitchDetector.h"
#include "CrossModHarmonicSeparator/HarmonicSeriesGenerator.h"
#include "CrossModHarmonicSeparator/HarmonicSpectralSeparator.h"


struct CrossModHarmonicSeparator
{
   static constexpr float kEnergyNormalizationGain = 1.0133f;

   // The UI elements defined in CrossModHarmonicSeparator.erbui are in 'ui'
   CrossModHarmonicSeparatorUi ui;

   // The data resources defined in CrossModHarmonicSeparator.erbb are in 'data'
   // Uncomment if you use some 'data' resources
   //CrossModHarmonicSeparatorData data;

   void  init ();
   void  process ();

   // Fundamental parameter setters (will map to UI controls in later phases)
   void setFundamentalLevel(float value);
   void setFundamentalTone(float value);
   void setFundamentalSpread(float value);

   // Prime parameter setters
   void setPrimeExciter(float value);
   void setPrimeBrightness(float value);
   void setPrimePhaseRotation(float value);
   // Composite parameter setters
   void setCompositeSaturation(float value);
   void setCompositeWarmth(float value);
   void setCompositeDensity(float value);

   const std::array<float, erb_BUFFER_SIZE>& fundamentalProcessedBlock() const
   {
      return fundamental_processed_block;
   }

   const std::array<float, erb_BUFFER_SIZE>& primeProcessedBlock() const
   {
      return prime_processed_block;
   }

   const std::array<float, erb_BUFFER_SIZE>& compositeProcessedBlock() const
   {
      return composite_processed_block;
   }

   // Put here your DSP objects
   PitchDetector             pitch_detector { 48000.0f };
   HarmonicSeriesGenerator   harmonic_series { 48000.0f };
   HarmonicSpectralSeparator spectral_separator { 48000.0f };
   FundamentalProcessor      fundamental_processor {};
   PrimeProcessor            prime_processor {};
   CompositeProcessor        composite_processor {};

   struct FundamentalParameters
   {
      float level = 1.0f;
      float tone = 0.5f;
      float spread = 0.0f;
   } fundamental_params;

   struct PrimeParameters
   {
      float exciter = 0.0f;
      float brightness = 0.5f;
      float phase_rotation = 0.0f;
   } prime_params;
   struct CompositeParameters
   {
      float saturation = 0.5f;
      float warmth = 0.5f;
      float density = 0.0f;
   } composite_params;

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
   std::array<float, erb_BUFFER_SIZE> fundamental_frame {};
   std::array<float, erb_BUFFER_SIZE> fundamental_processed_block {};
   std::array<float, erb_BUFFER_SIZE> prime_frame {};
   std::array<float, erb_BUFFER_SIZE> prime_processed_block {};
   std::array<float, erb_BUFFER_SIZE> composite_frame {};
   std::array<float, erb_BUFFER_SIZE> composite_processed_block {};
   
   // Debug members for energy analysis (only in test environment)
   #if defined(CROSSMOD_DEBUG_ENABLED) || defined(_DEBUG)
   std::vector<float> debug_window_sum_pattern;  // Sum of windows (not squared) for COLA verification
   std::vector<int> debug_overlap_count_pattern;  // Number of overlapping frames at each position
   double debug_pre_normalization_energy = 0.0;   // Energy before normalization
   double debug_post_normalization_energy = 0.0;   // Energy after normalization
   double debug_input_energy = 0.0;                // Input energy for current buffer
   std::vector<double> debug_normalization_values; // Normalization values used for output
   double debug_norm_min = 0.0;
   double debug_norm_max = 0.0;
   double debug_norm_avg = 0.0;
   double debug_spectral_ratio = 0.0;
   double debug_total_separated_energy = 0.0;
   double debug_energy_total_latest = 0.0;
   std::vector<float> debug_last_normalization_overlap;
   std::vector<float> debug_last_fundamental_overlap;
   std::vector<float> debug_last_prime_overlap;
   std::vector<float> debug_last_composite_overlap;
   std::vector<float> debug_last_input_buffer;
   void dumpDebugState(const char* label) const;
   #endif
};
