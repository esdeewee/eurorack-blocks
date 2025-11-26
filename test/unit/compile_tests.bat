cl /EHsc /std:c++17 /I. /I../../include /I../../submodules/vcv-rack-sdk/dep/include ^
main.cpp ^
../../src/detail/Animation.cpp ^
../../src/detail/Debounce.cpp ^
../../src/detail/Sdram.cpp ^
TestAnimation.cpp ^
TestDebounce.cpp ^
TestSdramPtr.cpp ^
frostflow/TestSpectralFft.cpp ^
../../samples/frostflow/dsp/SpectralFft.cpp ^
../../samples/frostflow/dsp/thirdparty/pffft.c ^
/Fe:test_runner.exe
