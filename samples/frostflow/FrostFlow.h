// FrostFlow.h

#pragma once

#include "artifacts/FrostFlowUi.h"
#include "erb/erb.h"
#include "FrostFlowDsp.h"

struct FrostFlow
{
   // The UI elements defined in FrostFlow.erbui are in 'ui'
   FrostFlowUi ui;
   
   // The DSP engine
   samples::frostflow::FrostFlowDsp dsp;

   void  init ();
   void  process ();
};
