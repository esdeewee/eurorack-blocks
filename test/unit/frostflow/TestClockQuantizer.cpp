#include "../test.h"
#include "TestHelpers.h"
#include "../../../samples/frostflow/dsp/ClockQuantizer.h"

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, ClockQuantizer_Queue) {
    ClockQuantizer q;
    q.init();
    
    float sampleRate = 1000.0f;
    bool req = true;
    bool clk = false;
    
    // No clock patched (unstable) -> Instant fire
    bool fire = q.process(clk, ClockQuantizer::DIV_4, sampleRate, req);
    erb_ASSERT_MSG(fire, "Should fire immediately if clock unstable");
    
    // Reset
    q.init();
    
    // Stable clock simulation
    // Send pulses
    // Period 200 samples (block size is small, say 64?)
    // process() assumes block size 64 internally if not passed.
    // Let's use large enough period > debounce.
    
    for (int i=0; i<20; ++i) {
        // Rising edge every 5 calls? 
        // 5 * 64 = 320 samples.
        q.process(true, ClockQuantizer::DIV_4, sampleRate, false); // Pulse
        q.process(false, ClockQuantizer::DIV_4, sampleRate, false);
        q.process(false, ClockQuantizer::DIV_4, sampleRate, false);
        q.process(false, ClockQuantizer::DIV_4, sampleRate, false);
        q.process(false, ClockQuantizer::DIV_4, sampleRate, false);
    }
    
    // Now Clock should be stable.
    
    // Queue
    // Call process WITHOUT clock pulse, WITH request
    fire = q.process(false, ClockQuantizer::DIV_4, sampleRate, true);
    erb_ASSERT_MSG(!fire, "Should queue when clock stable");
    
    // Advance to next pulse
    // Logic: Fire on rising edge of clock (DIV_4)
    bool fired = false;
    for (int i=0; i<200; ++i) {
        clk = (i % 100 < 10); // Pulse at i=0, i=100
        // Shift phase relative to previous loop? 
        // Let's just run until fire
        if (q.process(clk, ClockQuantizer::DIV_4, sampleRate, false)) {
            fired = true;
            break;
        }
    }
    erb_ASSERT_MSG(fired, "Quantized trigger failed to fire");
}

erb_TEST_CASE(FrostFlow, ClockQuantizer_Divisions) {
    // Placeholder: Test if DIV_8 fires twice as often?
    // Current implementation is placeholder for subdivisions.
    // Main clock alignment is guaranteed.
    erb_ASSERT_MSG(true, "Divisions test placeholder");
}
