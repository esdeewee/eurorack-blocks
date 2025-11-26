#include "../test.h"
#include "../../../samples/frostflow/dsp/SpectralProcessor.h"
#include <vector>

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, SpectralProcessor_Tilt) {
    SpectralProcessor proc;
    proc.init();
    
    std::vector<float> mag(512, 1.0f); 
    std::vector<float> phase(512, 0.0f);
    
    proc.setTilt(1.0f); 
    proc.process(mag.data(), phase.data(), 512);
    
    // Check low bin (DC) (should be attenuated)
    erb_TEST(mag[0] < 1.0f);
    
    // Check high bin (Nyquist) (should be boosted)
    erb_TEST(mag[1] > 1.0f);
}
