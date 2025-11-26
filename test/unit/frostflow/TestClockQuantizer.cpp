#include "../test.h"
#include "../../../samples/frostflow/dsp/ClockQuantizer.h"

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, ClockQuantizer_Queue) {
    ClockQuantizer q;
    q.init();
    
    float sampleRate = 1000.0f;
    
    // Queue a request
    bool request = true;
    bool clk = false;
    
    // First call queues it
    bool fire = q.process(clk, ClockQuantizer::DIV_4, sampleRate, request);
    
    // Next calls with no clock should fire immediately if unstable/no clock
    // My implementation: "if (!clockStable) ... if (triggerQueued) return true;"
    // Clock starts unstable.
    
    erb_TEST(fire == true); // Should fire immediately if no clock
}
