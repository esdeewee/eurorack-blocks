#include "../test.h"
#include "TestHelpers.h"
#include "../../../samples/frostflow/dsp/SidechainProcessor.h"

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, SidechainProcessor_Env) {
    SidechainProcessor sc;
    sc.init();
    
    bool trigger = false;
    float dryWet = 0.0f;
    
    // Attack phase
    // Need enough samples for envelope to rise
    // Attack coeff 0.1 -> reaches ~63% in 10 samples, ~99% in 45 samples.
    for(int i=0; i<100; ++i) {
        sc.process(1.0f, 1.0f, SidechainProcessor::ENV, dryWet, trigger);
    }
    
    float mod = sc.process(1.0f, 1.0f, SidechainProcessor::ENV, dryWet, trigger);
    
    erb_ASSERT_MSG(mod > 0.9f, "Envelope follower failed to rise fully");
    erb_ASSERT_MSG(mod <= 1.0f, "Envelope follower overshot");
}

erb_TEST_CASE(FrostFlow, SidechainProcessor_Gate) {
    SidechainProcessor sc;
    sc.init();
    
    bool trigger = false;
    float dryWet = 0.5f;
    
    // Low
    sc.process(0.0f, 1.0f, SidechainProcessor::GATE, dryWet, trigger);
    
    // High
    sc.process(1.0f, 1.0f, SidechainProcessor::GATE, dryWet, trigger);
    // Trigger should be true on rising edge
    erb_ASSERT_MSG(trigger, "Gate trigger failed");
    
    // Hold High
    sc.process(1.0f, 1.0f, SidechainProcessor::GATE, dryWet, trigger);
    erb_ASSERT_MSG(!trigger, "Gate re-triggered while holding");
}

erb_TEST_CASE(FrostFlow, SidechainProcessor_LFO) {
    SidechainProcessor sc;
    sc.init();
    bool trigger = false;
    float dryWet = 0.5f;
    
    // LFO Mode passes input CV directly
    // Input 1.0, Depth 0.5 -> Mod +0.5 -> 1.0
    float mod = sc.process(1.0f, 0.5f, SidechainProcessor::LFO, dryWet, trigger);
    erb_ASSERT_FLOAT_EQ(mod, 1.0f, 0.001f, "LFO mode positive scaling");
    
    // Input -1.0, Depth 0.5 -> Mod -0.5 -> 0.0
    mod = sc.process(-1.0f, 0.5f, SidechainProcessor::LFO, dryWet, trigger);
    erb_ASSERT_FLOAT_EQ(mod, 0.0f, 0.001f, "LFO mode negative scaling");
}
