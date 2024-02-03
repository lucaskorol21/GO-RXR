import os
import unittest

def load_tests_from_file(test_file):
    # Load the test suite from the given file
    loader = unittest.TestLoader()
    suite = loader.discover(start_dir=os.path.dirname(test_file), pattern=os.path.basename(test_file))
    return suite

def discover_and_run_tests(test_dir):
    # Discover and run all test files in the specified directory
    for root, _, files in os.walk(test_dir):
        for file in files:
            if file.startswith("test_") and file.endswith(".py"):
                test_file = os.path.join(root, file)
                print(f"Running tests from: {test_file}")
                test_suite = load_tests_from_file(test_file)
                runner = unittest.TextTestRunner()
                result = runner.run(test_suite)

    return result

if __name__ == "__main__":
    test_dir = "."
    result = discover_and_run_tests(test_dir)