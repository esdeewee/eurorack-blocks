#!/usr/bin/env bash
set -euo pipefail

MODULE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${MODULE_ROOT}/.." && pwd)"
BUILD_DIR="${MODULE_ROOT}/tests/build"
GTEST_ROOT="${REPO_ROOT}/submodules/libDaisy/tests/googletest"

mkdir -p "${BUILD_DIR}"

VCVARS_BAT="C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Auxiliary/Build/vcvars64.bat"

if [[ -f "${VCVARS_BAT}" ]]; then
    WIN_VCVARS=$(cygpath -w "${VCVARS_BAT}")
    WIN_BUILD_DIR=$(cygpath -w "${BUILD_DIR}")
    WIN_GTEST_INCLUDE=$(cygpath -w "${GTEST_ROOT}/googletest/include")
    WIN_REPO_ROOT=$(cygpath -w "${REPO_ROOT}")
    WIN_GTEST_ROOT=$(cygpath -w "${GTEST_ROOT}/googletest")
    WIN_HELPERS=$(cygpath -w "${MODULE_ROOT}/tests/helpers")
    WIN_MODULE_ROOT=$(cygpath -w "${MODULE_ROOT}")
    WIN_MODULE_CPP=$(cygpath -w "${MODULE_ROOT}/../CrossModHarmonicSeparator.cpp")
    WIN_HARMONIC_CPP=$(cygpath -w "${MODULE_ROOT}/HarmonicSeriesGenerator.cpp")
    WIN_SPECTRAL_CPP=$(cygpath -w "${MODULE_ROOT}/HarmonicSpectralSeparator.cpp")
    WIN_PITCH_CPP=$(cygpath -w "${MODULE_ROOT}/PitchDetector.cpp")
    WIN_KISSFFT_WRAPPER=$(cygpath -w "${MODULE_ROOT}/kissfft_wrapper.cpp")
    WIN_TEST_SIGNAL=$(cygpath -w "${MODULE_ROOT}/tests/test_signal_flow.cpp")
    WIN_TEST_PITCH=$(cygpath -w "${MODULE_ROOT}/tests/test_pitch_detection.cpp")
    WIN_TEST_HARMONIC=$(cygpath -w "${MODULE_ROOT}/tests/test_harmonic_series.cpp")
    WIN_GTEST_ALL=$(cygpath -w "${GTEST_ROOT}/googletest/src/gtest-all.cc")
    WIN_GTEST_MAIN=$(cygpath -w "${GTEST_ROOT}/googletest/src/gtest_main.cc")
    WIN_RESULTS_XML=$(cygpath -w "${BUILD_DIR}/test-results.xml")
    WIN_BATCH=$(cygpath -w "${BUILD_DIR}/run_tests_win.bat")

    cat > "${BUILD_DIR}/run_tests_win.bat" <<EOF
@echo off
call "${WIN_VCVARS}"
cl /nologo /std:c++17 /Zc:__cplusplus /EHsc ^
 /DCROSSMOD_TEST_ENV_HEADER=\"CrossModTestEnv.h\" ^
 /I"${WIN_GTEST_INCLUDE}" ^
 /I"${WIN_GTEST_ROOT}" ^
 /I"${WIN_REPO_ROOT}" ^
 /I"${WIN_HELPERS}" ^
 /I"${WIN_MODULE_ROOT}" ^
 "${WIN_MODULE_CPP}" ^
 "${WIN_HARMONIC_CPP}" ^
 "${WIN_SPECTRAL_CPP}" ^
 "${WIN_PITCH_CPP}" ^
"${WIN_KISSFFT_WRAPPER}" ^
 "${WIN_TEST_SIGNAL}" ^
 "${WIN_TEST_PITCH}" ^
 "${WIN_TEST_HARMONIC}" ^
 "${WIN_GTEST_ALL}" ^
 "${WIN_GTEST_MAIN}" ^
 /Fe:"${WIN_BUILD_DIR}\\runTests.exe"
if errorlevel 1 exit /b %errorlevel%
"${WIN_BUILD_DIR}\\runTests.exe" --gtest_output=xml:${WIN_RESULTS_XML}
exit /b %errorlevel%
EOF

    cmd.exe /C "\"${WIN_BATCH}\""
else
    CXX_COMPILER="${REPO_ROOT}/build-system/toolchain/msys2_mingw64/bin/g++.exe"
    if [[ ! -x "${CXX_COMPILER}" ]]; then
        CXX_COMPILER="g++"
    fi

    "${CXX_COMPILER}" -std=c++17 \
        -DCROSSMOD_TEST_ENV_HEADER=\"CrossModTestEnv.h\" \
        -I"${GTEST_ROOT}/googletest/include" \
        -I"${GTEST_ROOT}/googletest" \
        -I"${MODULE_ROOT}/.." \
        -I"${MODULE_ROOT}/tests/helpers" \
        -I"${MODULE_ROOT}" \
        "${MODULE_ROOT}/../CrossModHarmonicSeparator.cpp" \
        "${MODULE_ROOT}/HarmonicSeriesGenerator.cpp" \
        "${MODULE_ROOT}/HarmonicSpectralSeparator.cpp" \
        "${MODULE_ROOT}/PitchDetector.cpp" \
        "${MODULE_ROOT}/kissfft_wrapper.cpp" \
        "${MODULE_ROOT}/tests/test_signal_flow.cpp" \
        "${MODULE_ROOT}/tests/test_pitch_detection.cpp" \
        "${MODULE_ROOT}/tests/test_harmonic_series.cpp" \
        "${GTEST_ROOT}/googletest/src/gtest-all.cc" \
        "${GTEST_ROOT}/googletest/src/gtest_main.cc" \
        -o "${BUILD_DIR}/runTests"

    "${BUILD_DIR}/runTests" --gtest_output=xml:"${BUILD_DIR}/test-results.xml"
fi

echo "Tests completed. Results written to ${BUILD_DIR}/test-results.xml"
