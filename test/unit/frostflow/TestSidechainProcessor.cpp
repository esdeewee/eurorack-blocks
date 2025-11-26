#include "../test.h"
#include "../../../samples/frostflow/dsp/SidechainProcessor.h"

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, SidechainProcessor_Env) {
    SidechainProcessor sc;
    sc.init();
    
    bool trigger = false;
    // Env mode, depth 1.0
    float dryWet = 0.0f;
    // High input -> should modulate dryWet UP
    float mod = sc.process(1.0f, 1.0f, SidechainProcessor::ENV, dryWet, trigger);
    
    // modulation adds to dryWet. If input is 1.0, env follower rises.
    // after some samples
    for(int i=0; i<1000; ++i) {
        mod = sc.process(1.0f, 1.0f, SidechainProcessor::ENV, dryWet, trigger);
    }
    
    erb_TEST(mod > dryWet);
}
