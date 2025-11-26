<!-- ef26dc93-bd66-406e-9d8a-411f582832be a11ec5e1-4e66-4fb6-96b7-74aa116c0ca2 -->
# FrostFlow Deployment Plan

## Overview

Two-phase deployment: Software (VCV Rack simulator) first, then Hardware (Daisy firmware). Each phase includes comprehensive automated testing with zero manual intervention required.

---

## Phase 1: Software Development (VCV Rack Simulator)

### 1.1 Project Initialization

**Location:** `samples/frostflow/`

**Tasks:**

- Create project structure: `erbb init --name FrostFlow`
- Generate base files: `FrostFlow.erbb`, `FrostFlow.erbui`, `FrostFlow.h`, `FrostFlow.cpp`
- Configure build system: `erbb configure`
- Verify VCV Rack project generation in `artifacts/project_vcvrack.xcodeproj`

**Files to create:**

- `samples/frostflow/FrostFlow.erbb` - Build configuration
- `samples/frostflow/FrostFlow.erbui` - UI layout (12HP, kivu12 board, all controls)
- `samples/frostflow/FrostFlow.h` - Main module structure
- `samples/frostflow/FrostFlow.cpp` - Main processing loop
- `samples/frostflow/FrostFlowDsp.h` - DSP engine header
- `samples/frostflow/FrostFlowDsp.cpp` - DSP engine implementation

### 1.2 DSP Core Implementation

**1.2.1 FFT Wrapper Module**

**File:** `samples/frostflow/dsp/SpectralFft.h`, `SpectralFft.cpp`

**Implementation:**

- Wrap PFFFT library (from VCV Rack SDK) for 1024-point real FFT
- Implement overlap-add processing with Hann window
- Provide clean interface: `process(float* input, float* output, size_t size)`
- Handle FFT setup/teardown lifecycle

**Test Requirements:**

- Unit test: FFT round-trip accuracy (input → FFT → IFFT → output, verify <0.01% error)
- Unit test: Frequency response (sine wave at known frequency, verify correct bin)
- Unit test: Phase preservation (test with known phase relationships)
- Unit test: Windowing function (verify Hann window application)
- Performance test: FFT timing (<400μs per 1024-point transform)

**1.2.2 Spectral Freeze Engine**

**File:** `samples/frostflow/dsp/SpectralFreeze.h`, `SpectralFreeze.cpp`

**Implementation:**

- Two independent freeze buffers (Layer A, Layer B)
- Capture trigger logic (manual, gate, threshold, clock-quantized)
- Spectral buffer management (magnitude + simplified phase storage)
- Buffer state machine (idle, capturing, frozen, decaying)

**Test Requirements:**

- Unit test: Capture trigger conditions (verify all trigger modes work)
- Unit test: Buffer management (verify no buffer overflows, proper state transitions)
- Unit test: Threshold detection (test with signals above/below -40dB)
- Unit test: Clock quantization (verify timing accuracy, queue management)
- Integration test: Capture → freeze → playback cycle

**1.2.3 Blend System**

**File:** `samples/frostflow/dsp/BlendEngine.h`, `BlendEngine.cpp`

**Implementation:**

- Dry/Wet interpolation (0% = live, 100% = frozen)
- A/B interpolation (0% = Layer A, 100% = Layer B)
- Equal-power crossfade curves
- CV input processing and smoothing

**Test Requirements:**

- Unit test: Dry/Wet interpolation (test 0%, 25%, 50%, 75%, 100% positions)
- Unit test: A/B interpolation (test all blend positions)
- Unit test: Equal-power curves (verify power conservation)
- Unit test: Edge cases (0% dry/wet with A/B changes, 100% dry/wet)
- Unit test: CV scaling (verify CV input maps correctly to parameter range)

**1.2.4 Sidechain Processing**

**File:** `samples/frostflow/dsp/SidechainProcessor.h`, `SidechainProcessor.cpp`

**Implementation:**

- Envelope follower (for Env mode)
- Gate detector (for Gate mode)
- LFO extractor (for LFO mode)
- Depth scaling and modulation application

**Test Requirements:**

- Unit test: Envelope follower accuracy (test with known signals)
- Unit test: Gate detection (verify rising edge detection)
- Unit test: LFO extraction (test with low-frequency CV)
- Unit test: Depth scaling (verify bipolar depth control works correctly)

**1.2.5 Clock Quantization**

**File:** `samples/frostflow/dsp/ClockQuantizer.h`, `ClockQuantizer.cpp`

**Implementation:**

- Clock estimation (similar to `samples/kick/ClockEstimator.cpp`)
- Division calculations (1/4, 1/8, 1/16 notes)
- Queue management (one pending capture per layer)
- Timing accuracy

**Test Requirements:**

- Unit test: Clock estimation (test with various clock rates)
- Unit test: Division calculations (verify 1/4, 1/8, 1/16 math)
- Unit test: Queue management (test queue overflow, proper firing)
- Unit test: Timing accuracy (verify captures align to clock)

**1.2.7 Deterministic RNG (Crucial for Testing)**

**File:** `samples/frostflow/dsp/RngWrapper.h`

**Implementation:**
- Wrapper around standard RNG
- Must allow setting a specific seed
- Used by Motion LFO and any jitter/randomization

**1.2.8 Audio Safety**

**File:** `samples/frostflow/dsp/Limiter.h`

**Implementation:**
- Soft clip or lookahead limiter at output
- Prevents volume spikes from spectral accumulation


**File:** `samples/frostflow/dsp/SpectralProcessor.h`, `SpectralProcessor.cpp`

**Implementation:**

- Spectral tilt (frequency weighting for lows/highs)
- Decay envelope (infinite hold → quick fade)
- Motion LFO (random-walk animation of spectrum)

**Test Requirements:**

- Unit test: Spectral tilt (verify frequency weighting curves)
- Unit test: Decay curves (test infinite hold, various fade times)
- Unit test: Motion LFO (verify random-walk behavior, bounds checking)

### 1.3 UI Integration

**1.3.1 ERBUI Layout**

**File:** `samples/frostflow/FrostFlow.erbui`

**Controls to implement:**

- 9 Pots: Dry/Wet, A/B, Layer A Level, Layer B Level, Layer A Tilt, Layer B Tilt, Layer A Decay, Layer B Decay, SC Depth
- 2 Buttons: FREEZE A, FREEZE B
- 2 LEDs: FREEZE A LED, FREEZE B LED
- 3 Switches: SC MODE (3-position), DIV (3-position), MIRROR MODE (2-position)
- 5 CV Inputs: Dry/Wet CV, A/B CV, Layer A Level CV, Layer B Level CV, SC IN
- 1 CV Output: Mirror Out
- 3 Gate Inputs: FREEZE A Gate, FREEZE B Gate, CLOCK IN
- 2 Audio In, 2 Audio Out

**1.3.2 Main Module Structure**

**File:** `samples/frostflow/FrostFlow.h`, `FrostFlow.cpp`

**Implementation:**

- Module struct with `FrostFlowUi ui` and `FrostFlowDsp dsp`
- `init()` function for initialization
- `process()` function for audio processing
- Connect UI controls to DSP parameters
- Handle all capture modes
- Process sidechain modulation
- Generate mirror output

### 1.4 Automated Testing Suite (Software Phase)

**1.4.1 Unit Test Framework**

**Location:** `test/unit/frostflow/`

**Test Files:**

- `TestSpectralFft.cpp` - FFT accuracy and performance
- `TestSpectralFreeze.cpp` - Freeze capture and buffer management
- `TestBlendEngine.cpp` - Blend mathematics
- `TestSidechainProcessor.cpp` - Sidechain processing
- `TestClockQuantizer.cpp` - Clock quantization
- `TestSpectralProcessor.cpp` - Tilt, decay, motion

**Test Infrastructure:**

- Extend `test/unit/test.gyp` to include FrostFlow tests
- Use existing test framework (`test.h`, `erb_TEST` macros)
- All tests must be deterministic (no random seeds, fixed inputs)

**1.4.2 Integration Test Framework**

**Location:** `test/integration/frostflow/`

**Test Files:**

- `TestFullPipeline.cpp` - Complete audio pipeline (input → FFT → freeze → blend → IFFT → output)
- `TestStateMachine.cpp` - State transitions and mode switching
- `TestCvProcessing.cpp` - CV input/output processing

**Test Data:**

- Pre-recorded test signals: sine waves, white noise, real audio samples
- Test fixtures in `test/data/frostflow/`

**1.4.3 Performance Test Framework**

**Location:** `test/performance/frostflow/`

**Test Files:**

- `TestCpuLoad.cpp` - Measure CPU usage per buffer (must be <1ms)
- `TestMemoryUsage.cpp` - Verify SRAM usage (<128KB)
- `TestRealTimeGuarantees.cpp` - Worst-case scenario testing

**1.4.4 Audio Quality Test Framework**

**Location:** `test/audio/frostflow/`

**Test Files:**

- `TestFrequencyResponse.cpp` - Automated spectral analysis
- `TestArtifactDetection.cpp` - Aliasing, clicks, distortion detection
- `TestPhaseCoherence.cpp` - Phase relationship verification
- `TestGoldenMaster.cpp` - Bit-exact comparison against reference WAV with fixed RNG seed

**Test Tools:**

- Use FFT analysis for automated frequency response
- Compare input/output spectra
- Detect artifacts using spectral analysis
- Deterministic RNG seeding is mandatory for these tests

### 1.5 Build & Validation

**1.5.1 Build Process**

- `erbb configure` - Generate IDE projects
- `erbb build simulator` - Build VCV Rack plugin
- Verify build succeeds on all platforms (macOS, Windows, Linux)

**1.5.2 Automated Validation**

- Run all unit tests: `test/unit/run.py`
- Run integration tests
- Run performance tests (verify CPU <1ms, memory <128KB)
- Run audio quality tests (verify no artifacts)

**1.5.3 VCV Rack Testing**

- Load module in VCV Rack
- Verify all controls work
- Test audio processing
- Verify visual feedback (LEDs, etc.)

---

## Phase 2: Hardware Development (Daisy Firmware)

### 2.1 Hardware-Specific Adaptations

**2.1.1 FFT Library Selection**

- Use CMSIS DSP library (ARM-optimized) instead of PFFFT
- Implement FFT wrapper that works for both simulator (PFFFT) and hardware (CMSIS)
- Use conditional compilation: `#if defined(erb_TARGET_DAISY)`

**2.1.2 Memory Management**

**File:** `samples/frostflow/FrostFlowDsp.cpp`

**Implementation:**

- Use `erb::SramPtr` for spectral buffers (not SDRAM - fits in SRAM)
- Verify memory allocation in `init()`
- Use static assertions to verify module size <128KB

**2.1.3 Real-Time Constraints**

- Add CPU load monitoring using `daisy::CpuLoadMeter`
- Verify processing time <1ms per buffer
- Add safety checks for buffer underruns

### 2.2 Hardware Testing Suite

**2.2.1 Unit Tests (Hardware-Compatible)**

- Same unit tests as software phase
- Run on desktop (using UNIT_TEST macro)
- Verify algorithms work identically

**2.2.2 Hardware-Specific Tests**

**Location:** `test/hardware/frostflow/`

**Test Files:**

- `TestMemoryAllocation.cpp` - Verify SRAM allocation, no SDRAM needed
- `TestCpuLoadHardware.cpp` - Measure actual CPU usage on hardware
- `TestRealTimePerformance.cpp` - Verify no dropouts under load

**2.2.3 Firmware Build & Flash**

- `erbb build hardware` - Build Daisy firmware
- `erbb install` - Flash to hardware
- Verify firmware loads and runs

**2.2.4 Hardware Validation**

- Test all controls (pots, buttons, switches, CV inputs)
- Test audio I/O (verify signal levels, no clipping)
- Test CV outputs (verify Mirror Out voltage range)
- Test gate inputs (verify trigger detection)
- Test LEDs (verify visual feedback)

---

## Testing Strategy Details

### Automated Test Coverage Requirements

**Unit Tests:**

- 100% code coverage for all DSP algorithms
- All edge cases tested (0%, 50%, 100% positions, etc.)
- All state transitions tested
- All error conditions tested

**Integration Tests:**

- Full audio pipeline tested with various input signals
- State machine tested with all mode combinations
- CV processing tested with various CV ranges

**Performance Tests:**

- CPU usage measured and verified <1ms
- Memory usage measured and verified <128KB
- Real-time guarantees verified under worst-case conditions

**Audio Quality Tests:**

- Frequency response verified (automated spectral analysis)
- Artifacts detected automatically (aliasing, clicks, distortion)
- Phase coherence verified

### Test Execution Automation

**CI/CD Integration:**

- All tests run automatically on every commit
- Tests run on all platforms (macOS, Windows, Linux)
- Performance tests fail if CPU/memory limits exceeded
- Audio quality tests fail if artifacts detected

**Test Data Management:**

- Test signals stored in `test/data/frostflow/`
- Reference outputs stored for comparison
- Automated comparison with tolerance thresholds

---

## Success Criteria

### Phase 1 (Software) Complete When:

- All unit tests pass (100% coverage)
- All integration tests pass
- Performance tests verify CPU <1ms, memory <128KB
- Audio quality tests show no artifacts
- VCV Rack module loads and functions correctly
- All controls work as specified

### Phase 2 (Hardware) Complete When:

- All software tests still pass
- Hardware-specific tests pass
- Firmware builds and flashes successfully
- All hardware controls verified working
- Audio I/O verified (signal levels, no clipping)
- CV I/O verified (correct voltage ranges)
- Real-time performance verified (no dropouts)

---

## File Structure

```
samples/frostflow/
├── FrostFlow.erbb
├── FrostFlow.erbui
├── FrostFlow.h
├── FrostFlow.cpp
├── FrostFlowDsp.h
├── FrostFlowDsp.cpp
└── dsp/
    ├── SpectralFft.h
    ├── SpectralFft.cpp
    ├── SpectralFreeze.h
    ├── SpectralFreeze.cpp
    ├── BlendEngine.h
    ├── BlendEngine.cpp
    ├── SidechainProcessor.h
    ├── SidechainProcessor.cpp
    ├── ClockQuantizer.h
    ├── ClockQuantizer.cpp
    ├── SpectralProcessor.h
    └── SpectralProcessor.cpp

test/unit/frostflow/
├── TestSpectralFft.cpp
├── TestSpectralFreeze.cpp
├── TestBlendEngine.cpp
├── TestSidechainProcessor.cpp
├── TestClockQuantizer.cpp
└── TestSpectralProcessor.cpp

test/integration/frostflow/
├── TestFullPipeline.cpp
├── TestStateMachine.cpp
└── TestCvProcessing.cpp

test/performance/frostflow/
├── TestCpuLoad.cpp
├── TestMemoryUsage.cpp
└── TestRealTimeGuarantees.cpp

test/audio/frostflow/
├── TestFrequencyResponse.cpp
├── TestArtifactDetection.cpp
└── TestPhaseCoherence.cpp

test/hardware/frostflow/
├── TestMemoryAllocation.cpp
├── TestCpuLoadHardware.cpp
└── TestRealTimePerformance.cpp

test/data/frostflow/
├── sine_440hz.wav
├── white_noise.wav
└── [other test signals]
```

---

## Implementation Order

1. **Project Setup** - Initialize project structure
2. **FFT Wrapper** - Core spectral processing (with tests)
3. **Spectral Freeze** - Capture system (with tests)
4. **Blend Engine** - Mixing system (with tests)
5. **UI Integration** - Connect controls to DSP
6. **Sidechain** - Modulation system (with tests)
7. **Clock Quantization** - Timing system (with tests)
8. **Spectral Processing** - Tilt, decay, motion (with tests)
9. **Integration** - Full pipeline (with tests)
10. **Performance Optimization** - CPU/memory tuning
11. **Hardware Adaptation** - CMSIS DSP, memory management
12. **Hardware Testing** - Hardware-specific validation
13. **Final Validation** - Complete test suite execution

Each step includes comprehensive automated testing before proceeding to the next.

### To-dos

- [ ] Initialize FrostFlow project structure (erbb init, create base files)
- [ ] Implement FFT wrapper module with unit tests (SpectralFft.h/cpp)
- [ ] Implement spectral freeze engine with unit tests (SpectralFreeze.h/cpp)
- [ ] Implement blend system with unit tests (BlendEngine.h/cpp)
- [ ] Create ERBUI layout and connect UI controls to DSP
- [ ] Implement sidechain processing with unit tests (SidechainProcessor.h/cpp)
- [ ] Implement clock quantization with unit tests (ClockQuantizer.h/cpp)
- [ ] Implement spectral processing (tilt, decay, motion) with unit tests
- [ ] Create and run integration tests (full pipeline, state machine, CV processing)
- [ ] Create and run performance tests (CPU load, memory usage, real-time)
- [ ] Create and run audio quality tests (frequency response, artifacts, phase)
- [ ] Validate module in VCV Rack (all controls, audio processing, visual feedback)
- [ ] Adapt for hardware (CMSIS DSP, memory management, real-time constraints)
- [ ] Create and run hardware-specific tests (memory allocation, CPU load, real-time)
- [ ] Build and flash firmware to hardware, verify basic functionality
- [ ] Validate all hardware controls, audio I/O, CV I/O, LEDs, real-time performance
- [ ] Run complete test suite on hardware, verify all success criteria met