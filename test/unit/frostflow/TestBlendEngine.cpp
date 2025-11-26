#include "../test.h"
#include "../../../samples/frostflow/dsp/BlendEngine.h"

using namespace samples::frostflow;

erb_TEST_CASE(FrostFlow, BlendEngine_Basic) {
    BlendEngine engine;
    engine.init();
    
    float gainDry, gainA, gainB;
    
    // 1. Fully Dry
    engine.process(0.0f, 0.0f, gainDry, gainA, gainB);
    erb_TEST(std::abs(gainDry - 1.0f) < 0.001f);
    
    // 2. Fully Wet, A
    engine.process(1.0f, 0.0f, gainDry, gainA, gainB);
    erb_TEST(std::abs(gainDry - 0.0f) < 0.001f);
    erb_TEST(gainA > 0.5f); 
    erb_TEST(std::abs(gainB - 0.0f) < 0.001f);
    
    // 3. Fully Wet, B
    engine.process(1.0f, 1.0f, gainDry, gainA, gainB);
    erb_TEST(std::abs(gainA - 0.0f) < 0.001f);
    erb_TEST(gainB > 0.5f);
}
