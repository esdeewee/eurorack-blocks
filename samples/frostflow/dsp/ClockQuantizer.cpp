#include "ClockQuantizer.h"

namespace samples {
namespace frostflow {

void ClockQuantizer::init() {
    lastClockGate = false;
    samplesSinceLastClock = 0;
    clockPeriodSamples = 0;
    clockStable = false;
    triggerQueued = false;
}

bool ClockQuantizer::process(bool clockGate, Division division, float sampleRate, bool triggerRequest) {
    // Assume called once per block of ~32/64 samples.
    // We need to advance time by block size.
    // Let's define a fixed block size assumption or pass it.
    // For now, assuming 64 samples (erb_BUFFER_SIZE standard).
    // This is a rough approximation for control rate.
    const uint32_t BLOCK_SIZE = 64; 

    bool clockRisingEdge = clockGate && !lastClockGate;
    lastClockGate = clockGate;

    if (samplesSinceLastClock < 1000000) { // Prevent wrap-around safety
        samplesSinceLastClock += BLOCK_SIZE;
    }

    // Timeout (~2s)
    if (samplesSinceLastClock > (uint32_t)(2.0f * sampleRate)) {
        clockStable = false;
    }

    if (clockRisingEdge) {
        if (samplesSinceLastClock > 100 + BLOCK_SIZE) { // Debounce
            if (clockStable) {
                clockPeriodSamples = (clockPeriodSamples * 3 + samplesSinceLastClock) / 4;
            } else {
                clockPeriodSamples = samplesSinceLastClock;
                clockStable = true;
            }
        }
        samplesSinceLastClock = 0;
        
        // Fire on main clock edge if queued
        if (triggerQueued) {
            triggerQueued = false;
            return true;
        }
    }

    if (triggerRequest) {
        triggerQueued = true;
    }

    if (!clockStable) {
        if (triggerQueued) {
            triggerQueued = false;
            return true;
        }
        return false;
    }

    // Sub-division firing (Coarse block accuracy)
    // If we just passed a sub-division point in this block...
    // Simplified: Only support Main Clock alignment for now to ensure robustness.
    // Sub-divisions require precise phase tracking which is hard at block rate.
    
    return false;
}

}
}
