##############################################################################
#
#     test.gyp
#     Copyright (c) 2020 Raphael DINGE
#
#Tab=3########################################################################



{
   'conditions': [
      ['OS=="mac"', {
         'includes': [
            'xcode.gypi',
         ],
      }],

      ['OS=="linux"', {
         'includes': [
            'linux.gypi',
         ],
      }],
   ],

   'targets' : [
      {
         'target_name': 'test',
         'type': 'executable',

         'defines': [
            'erb_TARGET_UNIT_TEST',
            '_USE_MATH_DEFINES',
         ],

         'include_dirs': [
            '.',
            '../../include',
            '../../submodules/vcv-rack-sdk/dep/include',
            '../../samples/frostflow',
         ],

         'sources': [
            'main.cpp',
            'test.h',

            '../../include/erb/SdramObject.h',
            '../../include/erb/SdramObject.hpp',
            '../../include/erb/SdramPtr.h',
            '../../include/erb/SdramPtr.hpp',
            '../../include/erb/detail/Animation.h',
            '../../include/erb/detail/Animation.hpp',
            '../../include/erb/detail/Debounce.h',
            '../../include/erb/detail/MonotonicMemoryPool.h',
            '../../include/erb/detail/MonotonicMemoryPool.hpp',
            '../../include/erb/detail/Sdram.h',
            '../../include/erb/detail/Sdram.hpp',

            '../../src/detail/Animation.cpp',
            '../../src/detail/Debounce.cpp',
            '../../src/detail/Sdram.cpp',

            'TestAnimation.cpp',
            'TestAnimation.h',
            'TestDebounce.cpp',
            'TestDebounce.h',
            'TestSdramPtr.cpp',
            'TestSdramPtr.h',
            
            # FrostFlow DSP Sources
            '../../samples/frostflow/dsp/thirdparty/pffft.c',
            '../../samples/frostflow/dsp/SpectralFft.cpp',
            '../../samples/frostflow/dsp/SpectralFft.h',
            '../../samples/frostflow/dsp/SpectralFreeze.cpp',
            '../../samples/frostflow/dsp/SpectralFreeze.h',
            '../../samples/frostflow/dsp/BlendEngine.cpp',
            '../../samples/frostflow/dsp/BlendEngine.h',
            '../../samples/frostflow/dsp/SidechainProcessor.cpp',
            '../../samples/frostflow/dsp/SidechainProcessor.h',
            '../../samples/frostflow/dsp/ClockQuantizer.cpp',
            '../../samples/frostflow/dsp/ClockQuantizer.h',
            '../../samples/frostflow/dsp/SpectralProcessor.cpp',
            '../../samples/frostflow/dsp/SpectralProcessor.h',
            '../../samples/frostflow/FrostFlowDsp.cpp',

            # FrostFlow Tests
            'frostflow/TestSpectralFft.cpp',
            'frostflow/TestSpectralFreeze.cpp',
            'frostflow/TestBlendEngine.cpp',
            'frostflow/TestSidechainProcessor.cpp',
            'frostflow/TestClockQuantizer.cpp',
            'frostflow/TestSpectralProcessor.cpp',
            'frostflow/TestGoldenMaster.cpp',
            'frostflow/TestLimiter.cpp',
            'frostflow/TestIntegration.cpp',
            'frostflow/TestExtremeStress.cpp',
         ],

         'ldflags': [
            '-pthread',
         ],
      },
   ],
}
