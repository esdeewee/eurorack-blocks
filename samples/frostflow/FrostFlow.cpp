// FrostFlow.cpp

#include "FrostFlow.h"
#include <cstddef>

void  FrostFlow::init ()
{
   dsp.init();
}

void  FrostFlow::process ()
{
   // 1. Update UI Parameters
   
   dsp.paramDryWet = float(ui.dry_wet) + float(ui.dry_wet_cv);
   dsp.paramAbMorph = float(ui.ab_morph) + float(ui.ab_morph_cv);
   
   // Levels
   dsp.paramLayerALevel = float(ui.layer_a_level) + float(ui.layer_a_level_cv);
   dsp.paramLayerBLevel = float(ui.layer_b_level) + float(ui.layer_b_level_cv);
   
   // Tilt
   dsp.getLayerA().setTilt((float(ui.layer_a_tilt) - 0.5f) * 2.0f); // 0..1 -> -1..1
   dsp.getLayerB().setTilt((float(ui.layer_b_tilt) - 0.5f) * 2.0f);
   
   // Decay
   dsp.getLayerA().setDecay(float(ui.layer_a_decay));
   dsp.getLayerB().setDecay(float(ui.layer_b_decay));
   
   // Sidechain
   dsp.paramSidechainDepth = (float(ui.sc_depth) - 0.5f) * 2.0f;
   dsp.sidechainInput = float(ui.sc_in);

   if (ui.sc_mode.position_first()) {
       dsp.sidechainMode = samples::frostflow::SidechainProcessor::ENV;
   } else if (ui.sc_mode.position_center()) {
       dsp.sidechainMode = samples::frostflow::SidechainProcessor::GATE; 
   } else {
       dsp.sidechainMode = samples::frostflow::SidechainProcessor::LFO;
   }
   
   // Clock & Divisions
   dsp.clockGate = ui.clock_in.operator bool();
   
   if (ui.div_mode.position_first()) {
       dsp.clockDiv = samples::frostflow::ClockQuantizer::DIV_4;
   } else if (ui.div_mode.position_center()) {
       dsp.clockDiv = samples::frostflow::ClockQuantizer::DIV_8;
   } else {
       dsp.clockDiv = samples::frostflow::ClockQuantizer::DIV_16;
   }

   // Freeze Requests
   // dsp.freezeReqX should be true if button pressed or gate high
   dsp.freezeReqA = ui.freeze_a_gate.operator bool() || ui.freeze_a.held();
   dsp.freezeReqB = ui.freeze_b_gate.operator bool() || ui.freeze_b.held();

   // Update LEDs based on DSP internal freeze state
   // We need to expose isFrozenA/B from DSP to UI?
   // Currently they are static in process() which is hidden.
   // Let's move isFrozen state to class members in FrostFlowDsp.h so we can read them.
   // For now, just light up when button pressed (feedback is immediate).
   // Better: DSP should expose `bool isFrozenA` public member.
   // I will stick to button state for now to avoid header changes, but for production, DSP state feedback is better.
   
   // Actually, `dsp.freezeA` was my parameter, but now I use local static `isFrozenA`.
   // I should promote them to members.
   // I'll update FrostFlowDsp.h to make `bool isFrozenA` public readonly or getter.
   
   if (dsp.isFrozenA) ui.led_freeze_a.on(1.0f); else ui.led_freeze_a.off();
   if (dsp.isFrozenB) ui.led_freeze_b.on(1.0f); else ui.led_freeze_b.off();

   // 2. Process Audio
   dsp.process(&ui.in_l[0], &ui.in_r[0], &ui.out_l[0], &ui.out_r[0], erb_BUFFER_SIZE);
   
   // 3. Mirror Output
   // "MIRROR MODE switch (2-position: Dry/Wet/A/B)"
   // Switch has 2 positions? Code says 2.
   float mirrorVal = 0.0f;
   if (ui.mirror_mode.position_first()) {
       // Dry/Wet
       mirrorVal = dsp.paramDryWet;
   } else {
       // A/B
       mirrorVal = dsp.paramAbMorph;
   }
   // Output 0..1 -> -5V..+5V? 
   // CvOut Bipolar range is -1..1 (representing -5V..+5V usually).
   // Map 0..1 to -1..1
   ui.mirror_out = mirrorVal * 2.0f - 1.0f;
}
