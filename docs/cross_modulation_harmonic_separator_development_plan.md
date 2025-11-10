# Cross-Modulation Harmonic Separator - Development Plan

**Target Framework:** eurorack-blocks  
**Developer:** Cursor AI Agent  
**Approach:** Test-Driven Development with Evidence-Based Validation

---

## Module Overview

**Name:** Cross-Modulation Harmonic Separator  
**Type:** FX Module  
**Size:** 12HP recommended  
**Core Function:** Splits audio into three harmonic groups (Fundamental, Prime Harmonics, Composite Harmonics), processes each separately, and enables cross-modulation between groups.

**Key Features:**
- Fundamental path: level, tone control, spread (chorus effect)
- Prime harmonics path: exciter, brightness filter, phase rotation
- Composite harmonics path: saturation, warmth filter, diffusion density
- Cross-modulation matrix: envelope followers drive parameters across groups
- Display: 3x3 modulation matrix visualization

---

## Development Philosophy

1. **Build incrementally** - each component testable in isolation
2. **Test before integration** - verify each block before connecting
3. **Evidence-based validation** - every test produces measurable output
4. **No blind spots** - if it can't be tested automatically, redesign it

---

## Phase 1: Project Setup & Basic Signal Flow

### Step 1.1: Create Module Scaffold
**Task:** Initialize eurorack-blocks project structure

**Implementation:**
- Create `.erbui` file with basic module definition
- Define audio I/O: 1 stereo input, 1 stereo output
- Add placeholder controls (will populate later)
- Create `.erbb` build configuration
- Create `.cpp/.h` files for module logic
- Scaffold `tests/` directory with gtest harness, helper utilities, and `run_tests.sh`

**Testing Strategy:**
```
TEST: Project builds without errors
VALIDATION: Run `erbb configure` and `erbb build simulator`
EXPECTED: Zero compilation errors
EVIDENCE: Build log shows successful compilation
```

---

### Step 1.2: Implement Passthrough Audio
**Task:** Basic audio routing from input to output

**Implementation:**
- In `process()` function: copy input buffer to output buffer
- No processing yet, just signal flow verification

**Testing Strategy:**
```
TEST: Audio passthrough correctness
INPUT: 1kHz sine wave, -12dB
VALIDATION: 
  - Output matches input exactly (bit-perfect)
  - Calculate correlation coefficient (should be 1.0)
  - Measure RMS difference (should be < 0.001)
EXPECTED: Perfect passthrough
EVIDENCE: Correlation = 1.0, RMS diff ≈ 0
```

```
TEST: Latency measurement
INPUT: Impulse (single sample spike)
VALIDATION: Measure delay between input impulse and output impulse
EXPECTED: Latency = erb_BUFFER_SIZE samples (one buffer delay)
EVIDENCE: Exact sample count
```

---

## Phase 2: Harmonic Analysis & Separation

### Step 2.1: Fundamental Frequency Detection
**Task:** Implement pitch detection algorithm

**Implementation:**
- Use YIN (or McLeod) algorithm for F0 detection with deterministic parameters
- Store detected fundamental frequency per buffer
- Handle low-RMS buffers by returning last valid pitch (graceful silence handling)

**Testing Strategy:**
```
TEST: Known frequency detection accuracy
INPUT: Pure sine waves at known frequencies
  - 100Hz, 200Hz, 440Hz, 1000Hz, 2000Hz
VALIDATION:
  - Detected F0 within ±2Hz of actual frequency
  - Test at various amplitudes (-40dB to 0dB)
EXPECTED: >98% accuracy within ±2Hz
EVIDENCE: Log of detected vs actual frequencies
```

```
TEST: Harmonic signal detection
INPUT: Sawtooth wave (rich harmonics) at 200Hz
VALIDATION: Detected F0 = 200Hz ±3Hz
EXPECTED: Fundamental correctly identified despite harmonics
EVIDENCE: Detection log over 1000 buffers
```

```
TEST: Silence handling
INPUT: Zero signal
VALIDATION: No false detections, no crashes
EXPECTED: Detection returns "no pitch" state
EVIDENCE: State machine log shows proper null handling
```

---

### Step 2.2: Harmonic Series Generation
**Task:** Generate harmonic frequency series from detected fundamental

**Implementation:**
- Calculate harmonic frequencies: F0, 2F0, 3F0, 4F0, etc. up to Nyquist
- Categorize into three groups:
  - Fundamental: 1x F0
  - Primes: 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31... (up to audio range)
  - Composites: 4, 6, 8, 9, 10, 12, 14, 15, 16, 18, 20...

**Testing Strategy:**
```
TEST: Harmonic frequency calculation
INPUT: F0 = 100Hz
VALIDATION:
  - Harmonic 2 = 200Hz (prime)
  - Harmonic 3 = 300Hz (prime)
  - Harmonic 4 = 400Hz (composite)
  - Harmonic 5 = 500Hz (prime)
  - etc.
EXPECTED: Exact frequency multiples
EVIDENCE: Array of calculated frequencies matches expected
```

```
TEST: Prime number categorization
INPUT: Harmonic numbers 1-50
VALIDATION:
  - Primes correctly identified: 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47
  - Composites correctly identified: all others
  - Fundamental separate: 1
EXPECTED: 100% correct categorization
EVIDENCE: Boolean array matches mathematical prime test
```

```
TEST: Nyquist limit handling
INPUT: F0 = 10kHz (harmonics would exceed 24kHz Nyquist)
VALIDATION: 
  - No harmonics above 24kHz in list
  - Higher harmonics simply excluded, no errors
EXPECTED: Clean cutoff at Nyquist
EVIDENCE: Max frequency in array < 24000Hz
```

---

### Step 2.3: FFT-Based Harmonic Separation
**Task:** Use FFT analysis to isolate harmonic groups deterministically

**Implementation:**
- Perform a fixed-size FFT (e.g., 1024-point Hann window) on each buffer
- Identify the fundamental bin nearest to F0
- Build masks for prime-index harmonics and composite-index harmonics up to Nyquist
- Apply masks to the magnitude spectrum; optionally reconstruct time-domain signals via overlap-add using pre-allocated SDRAM buffers

**Testing Strategy:**
```
TEST: Fundamental bin identification
INPUT: Pure 440Hz sine
VALIDATION:
  - Detected peak bin matches theoretical index
  - Fundamental mask captures >95% of total spectral energy
EXPECTED: Accurate bin selection
EVIDENCE: Numerical comparison of bin magnitudes in unit test
```

```
TEST: Prime vs composite energy separation
INPUT: Sawtooth 200Hz
VALIDATION:
  - Prime mask contains energy at prime harmonic indices (2, 3, 5…)
  - Composite mask contains energy at composite indices (4, 6, 8…)
  - Crosstalk between masks <5%
EXPECTED: Clear spectral separation
EVIDENCE: Energy-per-mask report produced by automated test
```

```
TEST: Reconstruction error
INPUT: Complex signal (sum of sine, saw, noise)
VALIDATION:
  - Sum of masked inverse FFT outputs approximates original buffer
  - RMS error <10% (expected due to windowing)
EXPECTED: Controlled reconstruction quality
EVIDENCE: Integration test comparing original vs reconstructed buffers
```

---

### Phase 2 Manual Verification (VCV Rack)

> After the automated test suite is green, run this VCV checklist before unlocking Phase 3. It gives final confidence that the simulator matches expectations end-to-end.

**Before you start**
- Rebuild the plugin: `python build-system/scripts/erbb build simulator`.
- In VCV Rack, load a new empty patch and add the following helpers:
  - `Fundamental VCO` ×2 (or any basic oscillator)
  - `WT LFO` or similar as a quick modulation source (optional)
  - `VCA MIX` (4-channel mixer) or `Mix4`
  - `Scope` (for time-domain view)
  - `Spectra` or `Audible Instruments FFT` (for spectrum view)
  - `Audio-8` (audio interface to your speakers/headphones)
  - Our module `CrossModHarmonicSeparator`
- Connect `Audio-8` output 1→2 as usual so you can hear the result.
- Keep the Scope set to DC coupling, 5 ms/div, and normal amplitude gain to make “flat vs wiggly” differences obvious.

1. **Simulator bring-up**
   - Patch `VCO 1` sine into `CrossModHarmonicSeparator IN`.
   - Patch `CrossModHarmonicSeparator OUT` to a free channel on `VCA MIX`, then route that to `Audio-8` output 1 (so you can hear it) and to the `Scope`.
   - Patch the same output into the spectrum analyzer. Leave the analyzer on “Log amplitude” so harmonics are easy to see.
   - Confirm you hear a clean sine and no clicks/pops when toggling the module’s bypass (if you wired a bypass switch) or when reconnecting the cable.

2. **Pitch & harmonic sanity**
   1. **Single-tone check (440 Hz)**
      - Set VCO 1 frequency to 440 Hz (A4). In Fundamental VCO, that’s “NOTE A4” or “Frequency 440”.
      - Expectation: the Pitch readout (or debug log if you have it spitting to the console) stays within ±1 Hz of 440. On the spectrum analyzer you should see one tall spike at 440 Hz.
   2. **Rich-harmonic check (200 Hz saw)**
      - Change VCO 1 waveform to sawtooth at 200 Hz (set frequency knob to 200 Hz, waveform selector to saw).
      - Expectation: the detector still reports ~200 Hz. On the analyzer, primes (400 Hz, 600 Hz, …) should light up in the “prime” meter/debug output, composites (800 Hz, 1000 Hz, …) in the composite path, consistent with the automated energy tests.
   3. **Two-tone fundamental test (440 Hz + 660 Hz)**
      - Patch VCO 1 sine at 440 Hz into mixer channel 1.
      - Patch VCO 2 sine at 660 Hz into mixer channel 2.
      - Mix both signals together (set both channel gains to unity, output to our module input).
      - Expectation: Even though neither oscillator is running at 220 Hz, the detector should report ≈220 Hz. That happens because 440 and 660 share a common divisor (Greatest Common Divisor) of 220; it’s the true repeating period of the combined signal. On the analyzer, you’ll see peaks at 220 Hz, 440 Hz, 660 Hz, 880 Hz, etc. Our “fundamental” bucket should capture the 220 Hz bin, “prime” bucket should grab bins whose harmonic number is prime (2→440 Hz, 3→660 Hz, 5→1100 Hz…), and “composite” bucket the rest.

3. **Spectral separation (do de drie deel-signalen samen het origineel vormen?)**
   - *Waarom we dit doen:* onze code splitst het inkomende signaal in drie tijdreeksen (fundamental / primes / composites). Als we die drie weer bij elkaar optellen, moet je praktisch hetzelfde krijgen als het originele signaal. Dat bewijst dat de scheiding geen gekke gaten achterlaat.
   - *Wat je nodig hebt:* een mixer met minstens 4 kanalen (bv. `VCA MIX`), het `Scope`, en de drie debug-uitgangen `FUND`, `PRIME`, `COMP` die op het frontpaneel staan.
   - *Stap-voor-stap:*
     1. Patch de drie debug-uitgangen van `CrossModHarmonicSeparator` naar kanaal 1/2/3 van `VCA MIX`. Zet de kanaal- en mastergain van de mixer allemaal op 0 dB (unity).
     2. Patch de mixer-uitgang naar een **tweede** ingang van het Scope (bijvoorbeeld Scope Channel 2). Patch ook het originele ingangssignaal van de module naar Scope Channel 1.
     3. Zet het Scope op dubbele weergave (Channel 1 & 2 tegelijk). Channel 1 = origineel, Channel 2 = som van de drie deel-signalen. De golfvormen moeten vrijwel perfect over elkaar liggen; een klein verschil door venstering is oké.
     4. Laat één voor één (mute-knop of kabel eruit) een van de drie debug-kanalen wegvallen. Je hoort en ziet dan welk deel van het geluid ze bevatten – handig om te checken of fundamental/prime/composite logisch klinken.
   - *Reconstruction error check:* in de module hebben we een variabele `reconstruction_error` (wordt gelogd in de console wanneer je de simulator start). Tijdens deze test moet die onder 0.7 blijven. Zie je hogere waarden? Noteer het inkomende signaal en meld het voordat we verder gaan.
   - *Pop-test:* terwijl audio speelt, trek de kabel uit de module-ingang en steek hem terug, of toggle een eventuele bypass-schakelaar. Je mag geen klik of volume-sprong horen. Zo wel, dan eerst oplossen voordat we naar Phase 3 gaan.

4. **Edge-case sweep**
   - Silence: Turn both VCO levels down to zero or unplug them. All outputs of the module should drop to silence; spectrum should show only noise floor. Detector may hold the last pitch briefly but must not fabricate wild values.
   - Frequency sweep: Slowly turn VCO 1 frequency from ~50 Hz up to 2 kHz over ~10 seconds. Watch that detected pitch slides smoothly without stepping, and that prime/composite energy indicators don’t flicker erratically.
   - Noise injection: With VCO 1 at 200 Hz sine, use `WT LFO` (or `Noise` module) at very low level (≈−60 dB) into the mixer along with the sine. The prime/composite energy ratios should remain near their clean-signal values (within the ±0.02 tolerance enforced in our automated noise test).

5. **Sign-off**
   - Write down any anomalies (e.g., pitch jump, audible pop, misrouted harmonic energy).
   - Only proceed to Phase 3 when every bullet above behaves exactly as described.

---

## Phase 3: Per-Group Processing Chains

> **Library integration note:** Before starting Phase 3 implementation, add the MIT-licensed [SoundPipe](https://github.com/PaulBatchelor/Soundpipe) DSP library (e.g. as a git submodule or vendored dependency). We will rely on its filters, saturation, delays, and envelope utilities for the fundamental/prime/composite chains in subsequent phases (software simulator, tests, and Daisy hardware).

### Step 3.1: Fundamental Processing
**Task:** Implement fundamental path: level control, tone filter, spread effect

**Implementation:**
- Level: Simple gain multiplication
- Tone: Shelving EQ (boost/cut high frequencies)
- Spread: Chorus effect via short delay + LFO modulation with pre-allocated delay buffers (no runtime allocation)

**Testing Strategy:**
```
TEST: Level control linearity
INPUT: Sine 1kHz, -20dB
TEST POINTS: Level pot at 0%, 25%, 50%, 75%, 100%
VALIDATION:
  - Measure output RMS at each level
  - Verify linear relationship (dB scale)
EXPECTED: 
  - 0% = silence (−∞ dB)
  - 50% = -6dB from full
  - 100% = 0dB (unity gain)
EVIDENCE: Level vs output dB table
```

```
TEST: Tone control frequency response
INPUT: White noise
TEST POINTS: Tone at -100% (cut), 0% (flat), +100% (boost)
VALIDATION:
  - FFT comparison at each setting
  - -100%: High freq (>2kHz) attenuated >10dB
  - +100%: High freq (>2kHz) boosted >10dB
EXPECTED: Smooth shelving curve
EVIDENCE: Frequency response plots at 3 settings
```

```
TEST: Spread effect creates detuning
INPUT: Pure 440Hz sine
SETUP: Spread at 50%
VALIDATION:
  - FFT shows sidebands around 440Hz
  - Sidebands within ±10Hz of center
  - Creates "wider" sound (stereo width measurable)
EXPECTED: Chorus sidebands visible
EVIDENCE: FFT waterfall plot showing modulation sidebands
```

---

### Step 3.2: Prime Harmonics Processing
**Task:** Implement prime path: exciter, brightness filter, phase rotation

**Implementation:**
- Exciter: Even-harmonic waveshaper (`tanh`) with controllable gain
- Brightness filter: Second-order high-pass (`erb::Biquad`) with programmable cutoff
- Phase rotation: Cascade of precomputed all-pass sections (static coefficients, no runtime allocation)

**Testing Strategy:**
```
TEST: Exciter harmonic generation
INPUT: Pure 100Hz sine (fundamental)
SETUP: Exciter at 0%, 50%, 100%
VALIDATION:
  - FFT at 0%: Only 100Hz present
  - FFT at 50%: Harmonics appear but controlled
  - FFT at 100%: Strong harmonic content (200, 300, 500Hz since primes)
EXPECTED: Progressive harmonic enrichment
EVIDENCE: THD (Total Harmonic Distortion) measurement increases
```

```
TEST: Brightness filter slope
INPUT: White noise
SETUP: Brightness filter at various cutoff frequencies
VALIDATION:
  - Measure -3dB point
  - Measure rolloff slope (should be 12dB/octave or 24dB/octave)
  - Verify high frequencies pass, low frequencies attenuated
EXPECTED: Clean high-pass characteristic
EVIDENCE: Frequency response sweep plot
```

```
TEST: Phase rotation effect
INPUT: Complex signal (chord: 100Hz + 200Hz + 300Hz)
SETUP: Phase rotation at 0°, 90°, 180°, 270°
VALIDATION:
  - Measure phase relationships between components
  - At 0°: No change
  - At 90°: Each harmonic shifted 90° relative
  - Amplitude spectrum unchanged (only phase)
EXPECTED: Phase shifts without amplitude change
EVIDENCE: Phase spectrum analysis at each rotation setting
```

---

### Step 3.3: Composite Harmonics Processing
**Task:** Implement composite path: saturation, warmth filter, granular density

**Implementation:**
- Saturation: Soft clipping waveshaper
- Warmth filter: Low-pass filter with resonance
- Granular density: Replace true granular synthesis with a deterministic multi-tap diffuse delay network (fixed random seed) to simulate density while remaining CPU-friendly

**Testing Strategy:**
```
TEST: Saturation curve soft clipping
INPUT: Sine wave at increasing amplitudes (-40dB to +10dB)
SETUP: Saturation at 50%
VALIDATION:
  - Plot input vs output amplitude
  - At low levels: Near-linear
  - At high levels: Curve flattens (soft clip)
  - No hard clipping (no discontinuities)
EXPECTED: Smooth saturation curve, THD increases gradually
EVIDENCE: Transfer function plot + THD vs input level
```

```
TEST: Warmth filter response
INPUT: White noise
SETUP: Warmth filter cutoff swept 500Hz to 5kHz
VALIDATION:
  - Measure -3dB point at each setting
  - Low frequencies pass, high frequencies roll off
  - Resonance peak at cutoff (Q>1)
EXPECTED: Low-pass characteristic with warmth
EVIDENCE: Frequency response plot family
```

```
TEST: Diffuse delay density control
INPUT: Continuous 1kHz sine
SETUP: Density parameter from 0% to 100%
VALIDATION:
  - At 0%: Only direct signal (minimal diffusion)
  - At 100%: Multi-tap energy spread across delay taps (measured via correlation decrease)
  - RMS remains stable (±1dB)
EXPECTED: Controlled diffusion scaling
EVIDENCE: Automated test computing correlation and RMS for each setting
```

```
TEST: Diffuse delay smoothness
INPUT: Continuous signal with sudden density change
SETUP: Density jumps from 0% to 100%
VALIDATION:
  - Output remains click-free (no impulses detected)
  - Parameter smoothing keeps transitions <5ms ramp
EXPECTED: Smooth audible transition
EVIDENCE: Sample-by-sample check for discontinuities
```

---

## Phase 4: Envelope Followers

### Step 4.1: Implement Per-Group Envelope Detection
**Task:** Extract amplitude envelopes from each group's audio

**Implementation:**
- RMS or peak envelope follower per group
- Attack/release smoothing (fast attack, slower release)
- Output: 0.0 to 1.0 range per group

**Testing Strategy:**
```
TEST: Envelope follower response to impulse
INPUT: Single impulse (one-sample spike) in one group
VALIDATION:
  - Measure attack time (time to 90% of peak)
  - Measure release time (time from peak to 10%)
  - Expected: Attack < 5ms, Release 50-200ms
EXPECTED: Fast attack, smooth release
EVIDENCE: Envelope curve plot with time measurements
```

```
TEST: Envelope accuracy with known RMS
INPUT: Sine wave with known RMS values
TEST POINTS: -40dB, -20dB, -10dB, 0dB
VALIDATION:
  - Detected envelope matches theoretical RMS ±1dB
  - Linear relationship in dB domain
EXPECTED: Accurate RMS tracking
EVIDENCE: Table of input RMS vs detected envelope
```

```
TEST: Envelope independence between groups
INPUT: Signal only in Fundamental group (primes/composites silent)
VALIDATION:
  - Fundamental envelope responds
  - Prime envelope stays at 0.0
  - Composite envelope stays at 0.0
  - No crosstalk
EXPECTED: Perfect isolation
EVIDENCE: Envelope values logged over 1000 buffers
```

---

## Phase 5: Cross-Modulation Routing

### Step 5.1: Implement Modulation Matrix
**Task:** Route envelope followers to modulate other groups' parameters

**Fixed Routing:**
- Fundamental envelope → Primes phase rotation speed
- Prime envelope → Composites granular density
- Composite envelope → Fundamental spread amount

**Implementation:**
- Three modulation depth controls (one per routing)
- Envelope value (0-1) × depth → parameter modulation amount
- Apply first-order smoothing and hard clamps to keep parameters within safe musical ranges

**Testing Strategy:**
```
TEST: Fundamental → Prime modulation
INPUT: Amplitude burst in Fundamental group
SETUP: Modulation depth 100%
VALIDATION:
  - Prime phase rotation speed increases during burst
  - Modulation follows Fundamental envelope shape
  - Smoothly returns to baseline after burst
EXPECTED: Phase rotation tracks Fundamental amplitude
EVIDENCE: Plot showing Fundamental envelope + Prime phase rotation over time
```

```
TEST: Modulation depth scaling
INPUT: Steady signal in Fundamental
SETUP: Mod depth at 0%, 50%, 100%
VALIDATION:
  - At 0%: Prime phase rotation unaffected
  - At 50%: Half modulation range
  - At 100%: Full modulation range
  - Linear relationship
EXPECTED: Depth control scales modulation linearly
EVIDENCE: Modulation amount vs depth table
```

```
TEST: Prime → Composite modulation
INPUT: Varying signal in Prime group (0dB to silence to 0dB)
SETUP: Modulation depth 100%
VALIDATION:
  - Composite diffusion density changes with Prime level
  - At Prime silence: Density at minimum
  - At Prime loud: Density at maximum
EXPECTED: Density follows Prime envelope
EVIDENCE: Prime envelope + Composite diffusion parameter plot
```

```
TEST: Composite → Fundamental modulation
INPUT: Composite group envelope modulation
VALIDATION:
  - Fundamental spread amount tracks Composite envelope
  - Spread increases with Composite level
  - Returns to baseline when Composite silent
EXPECTED: Spread modulation works correctly
EVIDENCE: Composite envelope + Fundamental spread plot
```

```
TEST: Modulation smoothness (anti-zipper)
INPUT: Step change in modulating envelope (0.0 → 1.0 instantly)
VALIDATION:
  - Parameter change is smoothed over N samples
  - No audible clicks or zippers
  - Measure slew rate
EXPECTED: Smooth parameter interpolation
EVIDENCE: Parameter value plot showing smooth ramp
```

---

### Step 5.2: Test Cross-Modulation Interactions
**Task:** Verify complex multi-group interactions

**Testing Strategy:**
```
TEST: Cascading modulation
INPUT: Loud signal in Fundamental only
SETUP: All modulation depths 100%
VALIDATION:
  - Fundamental envelope high → Prime phase rotates fast
  - Prime processing creates strong Prime envelope
  - Prime envelope → Composite diffusion density increases
  - Composite processing creates Composite envelope
  - Composite envelope → Fundamental spread increases
  - Full chain verified
EXPECTED: Modulation cascade flows through all groups
EVIDENCE: Timeline plot showing all 3 envelopes + all modulated parameters
```

```
TEST: Feedback stability
INPUT: Strong signal cycling through groups with feedback routing
SETUP: Feedback loop: Fund → Prime → Comp → Fund (via spread)
VALIDATION:
  - System doesn't oscillate out of control
  - Parameters stay within valid ranges
  - No crashes or numerical overflow
EXPECTED: Stable feedback, controlled modulation
EVIDENCE: 10-second continuous run log showing parameter bounds
```

```
TEST: Modulation independence
INPUT: Silence in one group, signal in others
VALIDATION:
  - Modulation only occurs from active groups
  - Silent groups don't produce false modulation
EXPECTED: Conditional modulation (only when source active)
EVIDENCE: Modulation values logged with input states
```

---

## Phase 6: Mix & Output Stage

### Step 6.1: Implement Mix Controls
**Task:** Blend three groups back together with individual levels

**Implementation:**
- Three level controls (Fundamental, Primes, Composites)
- Crossfade between groups (optional)
- Master output level
- Dry/wet control (original signal vs processed)

**Testing Strategy:**
```
TEST: Individual group muting
INPUT: All groups active
TEST POINTS: Mute each group individually
VALIDATION:
  - Mute Fundamental: Only Primes + Composites audible
  - Mute Primes: Only Fund + Comp audible
  - Mute Composites: Only Fund + Prime audible
  - Mute all: Silence
EXPECTED: Complete isolation control
EVIDENCE: RMS measurement of output per test case
```

```
TEST: Mix level accuracy
INPUT: Known signal level in each group (-20dB)
SETUP: Mix levels at 0%, 50%, 100% for each group
VALIDATION:
  - Output level scales linearly with mix control
  - 0% = group absent
  - 100% = group at full level
EXPECTED: Linear level control
EVIDENCE: Output RMS vs mix level table per group
```

```
TEST: Dry/wet blend
INPUT: Complex signal
SETUP: Dry/wet at 0% (dry), 50%, 100% (wet)
VALIDATION:
  - 0%: Output = input (perfect passthrough)
  - 50%: Equal mix of dry and processed
  - 100%: Only processed signal
EXPECTED: Smooth crossfade, no level jumps
EVIDENCE: Correlation with input + RMS levels at each setting
```

```
TEST: Output clipping prevention
INPUT: All groups at maximum level simultaneously
VALIDATION:
  - Output never exceeds ±1.0
  - No hard clipping distortion
  - Soft limiting if needed
EXPECTED: Clean limiting, no distortion
EVIDENCE: Peak level measurement + THD check
```

---

## Phase 7: User Interface & Controls

### Step 7.1: Define Control Layout
**Task:** Map all controls to pots/encoders/buttons in .erbui

**Controls Needed (15 parameters):**
1. Input Gain
2. Fundamental: Level, Tone, Spread
3. Primes: Exciter, Brightness, Phase Rotation
4. Composites: Saturation, Warmth, Diffusion Density
5. Cross-Mod Depths: Fund→Prime, Prime→Comp, Comp→Fund
6. Mix: Fundamental, Primes, Composites
7. Dry/Wet

**Testing Strategy:**
```
TEST: Pot value reading
SETUP: Set each pot to 0%, 50%, 100% physically in simulator
VALIDATION:
  - Software reads correct values (±2% tolerance)
  - Full range 0.0 to 1.0 utilized
  - No dead zones at endpoints
EXPECTED: Accurate pot readings
EVIDENCE: Pot value log at each position
```

```
TEST: Control response latency
INPUT: Rapid pot change (0% to 100% in one buffer)
VALIDATION:
  - Parameter updates within erb_BUFFER_SIZE samples
  - No lag or sluggish response
EXPECTED: Immediate response
EVIDENCE: Timestamp log of pot change + parameter update
```

```
TEST: Encoder increments
SETUP: Rotary encoder for parameter selection
VALIDATION:
  - Clockwise: Increment by 1
  - Counter-clockwise: Decrement by 1
  - No missed steps
  - Boundaries respected (no overflow)
EXPECTED: Reliable encoder counting
EVIDENCE: Encoder count log during 100 rotations
```

---

### Step 7.2: Implement Display Visualization
**Task:** 128x64 OLED display showing modulation matrix

**Display Elements:**
- 3x3 grid representing modulation routing
- Bar graphs showing modulation depths
- Real-time envelope levels per group
- Active group indicators

**Testing Strategy:**
```
TEST: Display refresh rate
VALIDATION:
  - Display updates at reasonable rate (10-30 FPS)
  - No visible flicker
  - Smooth bar graph animations
EXPECTED: Fluid visual feedback
EVIDENCE: Frame timing log
```

```
TEST: Display data accuracy
SETUP: Set known modulation depths and envelope levels
VALIDATION:
  - Display shows correct values (±5% visual tolerance)
  - Envelope bars match actual envelope values
  - Modulation depth numbers match pot settings
EXPECTED: Accurate visual representation
EVIDENCE: Screenshot comparison with logged values
```

```
TEST: Display rendering performance
VALIDATION:
  - Display update doesn't cause audio glitches
  - Process timing stays within buffer deadline
  - Display runs in idle() thread, not process()
EXPECTED: Zero audio impact from display
EVIDENCE: Process timing log with/without display updates
```

---

### Step 7.3: Automated Test Harness Integration
**Task:** Wire up deterministic host-side tests for every DSP block

**Implementation:**
- Create `tests/` directory with gtest-based unit tests (mirroring `eurorack-blocks/test/unit`)
- Provide helper utilities (`generate_sine`, `generate_noise`, `compute_fft`, `compute_rms`) for repeatable DSP analysis
- Build tests as host executables via custom CMake or direct `g++` script invoked from Git Bash
- Supply `run_tests.sh` (cross-platform via Git Bash) that compiles and executes all tests, emitting JUnit XML for CI pipelines

**Testing Strategy:**
```
TEST: Test harness determinism
SETUP: Run entire test suite twice consecutively
VALIDATION:
  - All tests pass on both runs
  - Numerical outputs (JSON/CSV) identical between runs
EXPECTED: Fully deterministic behaviour
EVIDENCE: Script comparing run #1 vs run #2 artefacts
```

```
TEST: Test suite runtime budget
VALIDATION:
  - Full suite completes in <3 minutes on reference PC
  - Per-test timing recorded and reported
EXPECTED: CI-friendly execution time
EVIDENCE: Timing summary generated by `run_tests.sh`
```

```
TEST: Coverage instrumentation
VALIDATION:
  - Tests built with coverage flags (gcov or llvm-cov)
  - Coverage report demonstrates >90% statement coverage for module sources
EXPECTED: High coverage achieved
EVIDENCE: Coverage summary stored in `artifacts/coverage/`
```

---

## Phase 8: Integration Testing

### Step 8.1: Full Signal Chain Test
**Task:** Verify complete audio path with all processing active

**Testing Strategy:**
```
TEST: End-to-end processing
INPUT: Recorded musical material (drums, guitar, synth)
SETUP: All groups active, moderate settings
VALIDATION:
  - Output is musically coherent
  - No glitches, clicks, pops
  - No unexpected silence or dropouts
  - Modulation creates audible changes
EXPECTED: Professional quality processing
EVIDENCE: Integration test computing THD, spectral centroid, and group energy vs expected thresholds
```

```
TEST: Extreme settings stability
SETUP: All parameters at maximum simultaneously
VALIDATION:
  - No crashes
  - No numerical overflows (NaN, Inf)
  - Output stays within ±1.0
  - Audio may be extreme but stable
EXPECTED: Robust handling of extreme values
EVIDENCE: 60-second simulated render with assertions on finite sample counts and bounded output
```

```
TEST: Silence handling
INPUT: Zero signal (complete silence)
VALIDATION:
  - Output is silent (not noise)
  - No self-oscillation
  - All envelopes at zero
  - Display shows zero activity
EXPECTED: Clean silence handling
EVIDENCE: Automated check confirming RMS < -100dB and all envelopes <0.01
```

---

### Step 8.2: Performance & Resource Testing

**Testing Strategy:**
```
TEST: CPU usage measurement
VALIDATION:
  - Measure process() execution time per buffer
  - Must complete within buffer deadline (erb_BUFFER_SIZE / 48000 seconds)
  - Typically should use <80% of available time
EXPECTED: Real-time performance maintained
EVIDENCE: Instrumented performance log parsed automatically to confirm average <80% and peak <95%
```

```
TEST: Memory allocation check
VALIDATION:
  - All large buffers (delays, FFT, etc.) allocated once in init()
  - No allocations in process() or idle()
  - Memory usage stable (no leaks)
  - SDRAM usage logged
EXPECTED: Static memory footprint
EVIDENCE: Unit test verifying allocator counters unchanged across 10,000 buffers
```

```
TEST: Simulator vs Hardware consistency
VALIDATION:
  - Feed identical deterministic test vectors to simulator and hardware
  - Collect numeric metrics (RMS, spectral centroid, latency) from both platforms
  - Differences stay within predefined tolerance bands
EXPECTED: Cross-platform consistency
EVIDENCE: Script comparing simulator metric JSON vs hardware UART metric log within tolerance
```

---

## Phase 9: Edge Cases & Robustness

### Step 9.1: Boundary Condition Testing

**Testing Strategy:**
```
TEST: Very low input levels
INPUT: -80dB sine wave
VALIDATION:
  - Pitch detection fails gracefully (no false positives)
  - Processing doesn't introduce noise floor
  - Envelopes return to zero properly
EXPECTED: Clean low-level handling
EVIDENCE: Noise floor measurement <-90dB
```

```
TEST: Very high input levels (clipping input)
INPUT: Clipped waveform (peaks at ±1.0)
VALIDATION:
  - Pitch detection still works or fails gracefully
  - No crashes from unexpected waveform shapes
  - Output limited properly
EXPECTED: Robust overload handling
EVIDENCE: Pitch detection log, output peak levels
```

```
TEST: Rapid parameter changes
SETUP: Deterministic pseudo-random sequence (fixed seed) updates all pot values every buffer
VALIDATION:
  - No audio glitches from rapid changes
  - Smoothing prevents zippering
  - System remains stable
EXPECTED: Smooth parameter interpolation handles chaos
EVIDENCE: Automated torture test logging max delta of output samples and confirming finite values
```

```
TEST: Extreme pitch detection cases
INPUT: Very low frequencies (20Hz) and very high (10kHz)
VALIDATION:
  - 20Hz: Detection works or fails gracefully
  - 10kHz: Detection works or fails gracefully
  - No crashes at extremes
EXPECTED: Robust across full audio range
EVIDENCE: Detection success rate log per frequency band
```

---

## Phase 10: Documentation & Finalization

### Step 10.1: Create Automated Test Suite
**Task:** Package all tests into executable test suite

**Implementation:**
- Provide `run_tests.sh` (Git Bash) that builds/executes host gtest binaries and parses performance logs
- Ensure script outputs JUnit XML, coverage summaries, and metric JSON for CI consumption
- Integrate script into automated pipeline (e.g., GitHub Actions) invoking Git Bash on Windows

**Testing Strategy:**
```
TEST: Test suite completeness
VALIDATION:
  - All module features have ≥1 test
  - Test coverage report generated
  - All critical paths tested
EXPECTED: >90% code coverage
EVIDENCE: Coverage report generated by `run_tests.sh`
```

```
TEST: Test suite reliability
VALIDATION:
  - Run test suite 10 times
  - All tests pass every time (deterministic)
  - No flaky tests
  - Total run time <3 minutes
EXPECTED: 100% pass rate, repeatable
EVIDENCE: Ten consecutive run logs archived by CI pipeline
```

---

### Step 10.2: User Documentation
**Task:** Create user manual and technical documentation

**Contents:**
- Module overview and concept
- Control descriptions
- Patch ideas and usage examples
- Technical specifications
- Troubleshooting guide

**No automated testing needed** (documentation quality review is manual)

---

### Step 10.3: Build & Flash to Hardware
**Task:** Final build for Daisy hardware

**Implementation:**
- Configure for `section qspi` if needed (check program size)
- Build firmware: `erbb build`
- Install bootloader (if needed): `erbb install bootloader`
- Flash to hardware: `erbb install firmware`

**Testing Strategy:**
```
TEST: Firmware size check
VALIDATION:
  - Program fits in target memory section
  - FLASH usage <100% (or QSPI <100% if using)
  - All resources allocated successfully
EXPECTED: Successful linking with memory headroom
EVIDENCE: Build log memory report
```

```
TEST: Hardware boot test
VALIDATION:
  - Module boots on Daisy hardware
  - Display initializes and shows content
  - Audio passes through (basic test)
  - No crashes on startup
EXPECTED: Clean boot sequence
EVIDENCE: Serial log from hardware boot
```

```
TEST: Hardware audio quality
INPUT: Test signals through hardware module
VALIDATION:
  - Match simulator output (within tolerance)
  - No unexpected noise or artifacts
  - Same latency as simulator
EXPECTED: Hardware = simulator quality
EVIDENCE: Automated metric comparison between simulator CSV and hardware UART output (RMS, spectral centroid, latency)
```

---

## Development Order Summary

**Recommended Build Sequence:**

1. **Foundation** (Phase 1): Project setup → Passthrough audio
2. **Analysis** (Phase 2): Pitch detection → Harmonic generation → Separation
3. **Processing** (Phase 3): Per-group chains (Fundamental → Primes → Composites)
4. **Modulation** (Phase 4-5): Envelope followers → Cross-mod routing
5. **Output** (Phase 6): Mix stage → Dry/wet
6. **Interface** (Phase 7): Controls → Display
7. **Integration** (Phase 8): Full chain testing → Performance
8. **Polish** (Phase 9-10): Edge cases → Documentation → Hardware

**Critical Path Items:**
- Harmonic separation quality determines everything else
- Envelope followers must be accurate for good modulation
- CPU performance must be monitored continuously (avoid late optimization)

---

## Testing Metrics Summary

### Automated Test Coverage

| Component | Test Count | Coverage Target |
|-----------|------------|-----------------|
| Signal Flow | 2 | 100% |
| Pitch Detection | 4 | ±2 Hz accuracy |
| Harmonic Separation | 3 | >95% energy isolation |
| Per-Group Processing | 6 | >90% |
| Envelope Followers | 3 | 100% |
| Cross-Modulation | 4 | 100% |
| Mix/Output | 3 | 100% |
| UI & Display | 4 | 100% |
| Integration | 3 | Full system |
| Edge Cases | 3 | Robustness |
| **TOTAL** | **35 automated tests** | **>90% coverage** |

### Test Evidence Types

- **Numerical**: RMS, peak levels, frequencies, timing
- **Comparative**: Input vs output correlation, difference
- **Boolean**: Pass/fail conditions, state checks
- **Visual**: Plots, spectrograms, frequency responses (auto-generated)
- **Statistical**: Distribution analysis, probability verification

### Test Automation Requirements

- All tests must be **executable via command line**
- All tests must produce **machine-readable output** (JSON/CSV)
- All tests must have **clear pass/fail criteria** (no human judgment)
- Test suite must run **in CI/CD environment** via Git Bash invoking host gtest binaries (simulator/performance builds optional but scriptable)

---

## Development Tools

### Required Tools
- **eurorack-blocks** framework
- **VCV Rack** simulator for testing
- **Test audio files** (sine waves, noise, musical content)
- **Python/C++ test harness** for automated validation
- **gtest** (GoogleTest) for C++ unit testing
- **Git Bash** (on Windows) to run build/test scripts
- **FFT analysis library** for spectral tests
- **SoundPipe** DSP library (MIT) for reusable building blocks in Phase 3+
- **Plotting library** for visualization (matplotlib, gnuplot)

### Test Data Generation
Automated tests synthesize signals in C++ for determinism. The following optional Python snippets can be used for manual validation or external tooling.
```python
# Example: Generate test signals
import numpy as np
import soundfile as sf

# Pure sine at 440Hz
t = np.linspace(0, 1, 48000)
sine_440 = np.sin(2 * np.pi * 440 * t)
sf.write('test_sine_440hz.wav', sine_440, 48000)

# Sawtooth with harmonics
sawtooth = 2 * (t % (1/100)) * 100 - 1
sf.write('test_sawtooth_100hz.wav', sawtooth, 48000)

# White noise
noise = np.random.normal(0, 0.1, 48000)
sf.write('test_white_noise.wav', noise, 48000)
```

### Test Validation Scripts
```python
# Example: Validate pitch detection
import json

with open('pitch_detection_results.json', 'r') as f:
    results = json.load(f)

for test in results['tests']:
    expected = test['expected_freq']
    detected = test['detected_freq']
    error = abs(detected - expected)
    
    assert error < 2.0, f"Pitch detection error {error}Hz exceeds tolerance"
    
print("All pitch detection tests PASSED")
```

---

## Risk Mitigation

### Technical Risks

| Risk | Mitigation Strategy | Test Coverage |
|------|---------------------|---------------|
| Pitch detection failures | Multiple algorithms, graceful degradation | 4 tests |
| CPU overload | Profile early, optimize hot paths | Continuous monitoring |
| Cross-mod instability | Bounded parameters, anti-runaway limiters | Stress tests |
| Harmonic crosstalk | High-Q filters, test isolation rigorously | 3 tests |
| Memory fragmentation | Static allocation only, no runtime malloc | Memory audit |

---

## Success Criteria

**Module is complete when:**

1. ✅ All 35 automated tests pass
2. ✅ CPU usage <80% in worst case
3. ✅ Memory allocation static and within limits
4. ✅ Hardware and simulator produce identical output
5. ✅ Musical testing confirms expected behavior
6. ✅ No crashes in 1-hour stress test
7. ✅ CI pipeline runs `run_tests.sh` with green results
8. ✅ Documentation complete

**Module is excellent when:**

1. ⭐ Modulation creates obvious and interesting sonic changes
2. ⭐ Cross-group interactions produce unexpected but musical results
3. ⭐ All controls have meaningful audible impact
4. ⭐ Display visualization is informative and beautiful
5. ⭐ Performance has headroom for additional features

---

## Cursor AI Specific Instructions

### Development Workflow

1. **Start with tests**: Write test before implementing feature
2. **Implement minimally**: Make test pass with simplest code
3. **Verify evidence**: Check test output matches expected evidence
4. **Refactor if needed**: Improve code while keeping tests green
5. **Run automation**: Execute `run_tests.sh` before promoting a change
6. **Document**: Add comments explaining non-obvious logic

### Code Organization

```
module/
├── CrossModHarmonicSeparator.erbui     # UI definition
├── CrossModHarmonicSeparator.erbb      # Build config
├── CrossModHarmonicSeparator.cpp       # Main module logic
├── CrossModHarmonicSeparator.h         # Module header
├── PitchDetector.h                     # Pitch detection algorithm
├── HarmonicSeparator.h                 # FFT-based spectral masking
├── EnvelopeFollower.h                  # Amplitude tracking
├── CrossModMatrix.h                    # Modulation routing
├── tests/
│   ├── test_signal_flow.cpp
│   ├── test_pitch_detection.cpp
│   ├── test_harmonic_separation.cpp
│   ├── test_processing_chains.cpp
│   ├── test_cross_modulation.cpp
│   ├── test_mix_output.cpp
│   └── test_integration.cpp
└── tests/helpers/
    ├── DSPTestUtils.h
    └── TestSignalGenerators.h
```

### Implementation Priorities

**Phase 1-2 (Foundation)**: Focus on correctness over performance
**Phase 3-5 (Processing)**: Balance correctness and performance
**Phase 6-7 (Polish)**: Optimize hot paths, refine UX
**Phase 8-10 (Completion)**: Robustness and documentation

### When Stuck

1. **Re-read test**: What exactly is being validated?
2. **Isolate component**: Test in isolation before integration
3. **Check examples**: Look at existing eurorack-blocks samples
4. **Simplify**: Remove complexity until test passes, then add back
5. **Verify assumptions**: Print intermediate values, check math

---

## Conclusion

This development plan provides a **complete roadmap** for implementing the Cross-Modulation Harmonic Separator module using **test-driven development** principles. Every feature has **automated validation**, ensuring **zero blind spots** in the implementation.

The plan is **AI-friendly**: each step is **clearly defined**, **independently testable**, and builds **incrementally** toward the final module. Cursor AI can follow this plan **phase by phase**, verifying each component before moving forward.

**Expected Development Time**: 20-40 hours of focused development
**Expected Code**: ~2000-3000 lines C++ (including tests)
**Expected Test Suite Run Time**: <3 minutes for full validation

**The result**: A robust, well-tested, musically interesting FX module ready for commercial production.
