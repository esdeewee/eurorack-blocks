import os
import subprocess
import sys
import re

def get_test_list():
    tests = []
    test_dir = "frostflow"
    if not os.path.exists(test_dir):
        print(f"Error: {test_dir} not found. Run from test/unit/")
        return []
        
    for f in os.listdir(test_dir):
        if f.endswith(".cpp"):
            with open(os.path.join(test_dir, f), 'r') as src:
                content = src.read()
                matches = re.findall(r'erb_TEST_CASE\s*\(\s*(\w+)\s*,\s*(\w+)\s*\)', content)
                for m in matches:
                    tests.append(f"{m[0]}::{m[1]}")
    return tests

def run_tests():
    # Path to test executable (Adjust based on actual build output)
    # MSBuild output: test/unit/Default/test.exe
    exe_path = os.path.abspath(os.path.join("Default", "test.exe"))
    
    if not os.path.exists(exe_path):
        print(f"Executable not found at {exe_path}")
        print("Please build first: python ../../build-system/scripts/erbb build unit-test")
        return

    tests = get_test_list()
    print(f"Found {len(tests)} tests.")
    
    passed = 0
    failed = []
    
    for test_name in tests:
        print(f"[EXEC] {test_name:<40} ", end='', flush=True)
        
        try:
            # Run the test executable with the specific test filter
            result = subprocess.run(
                [exe_path, test_name],
                capture_output=True,
                text=True,
                timeout=10 # Increased timeout for torture tests
            )
            
            if result.returncode == 0:
                print("PASS")
                passed += 1
            else:
                print(f"FAIL (Exit Code: {result.returncode})")
                print("--- Output ---")
                print(result.stdout)
                print("--- Error ---")
                print(result.stderr)
                failed.append(test_name)
                
        except subprocess.TimeoutExpired:
            print("TIMEOUT")
            failed.append(test_name)
        except Exception as e:
            print(f"CRASH: {e}")
            failed.append(test_name)

    print("-" * 60)
    print(f"Summary: {passed}/{len(tests)} Passed")
    if failed:
        print("Failed Tests:")
        for t in failed:
            print(f" - {t}")
        sys.exit(1)
    else:
        print("All tests passed successfully.")

if __name__ == "__main__":
    run_tests()
