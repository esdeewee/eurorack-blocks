#pragma once

#include <cstdint>

namespace samples {
namespace frostflow {

class ClockQuantizer {
public:
    enum Division {
        DIV_4,  // Quarter note (1 beat)
        DIV_8,  // Eighth note
        DIV_16  // Sixteenth note
    };

    void init();

    /**
     * Process clock input and trigger queue.
     * @param clockGate Gate input state of the clock
     * @param division Clock division setting
     * @param sampleRate Audio sample rate
     * @param triggerRequest Input request to trigger a freeze (e.g. from button/gate)
     * @return True if a quantized trigger should fire now
     */
    bool process(bool clockGate, Division division, float sampleRate, bool triggerRequest);

private:
    // Clock Tracking
    bool lastClockGate = false;
    uint32_t samplesSinceLastClock = 0;
    uint32_t clockPeriodSamples = 0; // Measured period
    bool clockStable = false; // Do we have a valid clock?

    // Division Logic
    uint32_t subClockCounter = 0;
    
    // Queue
    bool triggerQueued = false;
};

}
}
