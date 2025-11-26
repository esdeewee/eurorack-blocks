#pragma once

#include <vector>
#include <string>
#include <functional>
#include <iostream>
#include <exception>
#include <cstdio>

// --- Assertion Macros ---
#if defined (__GNUC__) || defined (__clang__)
   #define erb_ASSERT(expr) \
      do { if (!(expr)) {printf ("      FAILED: %s line: %d\n      Expression: %s\n", __PRETTY_FUNCTION__, __LINE__, #expr); std::terminate ();}} while (false)
#elif defined (_MSC_VER)
   #define erb_ASSERT(expr) \
      do { if (!(expr)) {printf ("      FAILED: %s line: %d\n      Expression: %s\n", __FUNCDNAME__, __LINE__, #expr); fflush(stdout); std::terminate ();}} while (false)
#else
   #error Unsupported Compiler.
#endif

// Alias for backward compatibility if needed, though I used erb_ASSERT in my new files.
// The original file had erb_TEST as the assertion.
#define erb_TEST(expr) erb_ASSERT(expr)


// --- Test Registry System ---

struct TestCase {
    std::string group;
    std::string name;
    std::function<void()> func;
};

class TestRegistry {
public:
    static TestRegistry& instance() {
        static TestRegistry inst;
        return inst;
    }

    void registerTest(const std::string& group, const std::string& name, std::function<void()> func) {
        tests.push_back({group, name, func});
    }

    void runAll() {
        bool allPassed = true;
        for (const auto& test : tests) {
            std::cout << "[RUNNING] " << test.group << "::" << test.name << " ... " << std::flush;
            try {
                test.func();
                std::cout << "[PASS]" << std::endl;
            } catch (...) {
                std::cout << "[FAIL] Exception caught" << std::endl;
                allPassed = false;
                // We continue? std::terminate in assert kills process.
                // To support crash-safe runner, we want the process to die on assertion failure?
                // Yes, the runner script handles the restart.
                // So we don't need to catch here if we rely on runner.
            }
        }
    }
    
    // Run matching filter
    void runFiltered(const std::string& filter) {
        for (const auto& test : tests) {
            std::string fullName = test.group + "::" + test.name;
            if (filter.empty() || fullName.find(filter) != std::string::npos) {
                std::cout << "[RUNNING] " << fullName << " ... " << std::flush;
                test.func();
                std::cout << "[PASS]" << std::endl;
            }
        }
    }

private:
    std::vector<TestCase> tests;
};

struct AutoRegister {
    AutoRegister(const std::string& group, const std::string& name, std::function<void()> func) {
        TestRegistry::instance().registerTest(group, name, func);
    }
};

#define erb_TEST_CASE(group, name) \
    static void group##_##name(); \
    static AutoRegister reg_##group##_##name(#group, #name, group##_##name); \
    static void group##_##name()
