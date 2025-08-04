#!/usr/bin/env python3

"""
Comprehensive test runner for the gene curation pipeline.

Runs all unit tests and generates a coverage report. Supports various
output formats and filtering options.
"""

import unittest
import sys
import os
import argparse
import time
from typing import List, Optional

# Add current directory to path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

def discover_tests(test_dir: str = "gene_curation_pipeline/tests",
                  pattern: str = "test_*.py") -> unittest.TestSuite:
    """Discover all test modules."""
    loader = unittest.TestLoader()
    
    # Convert relative path to absolute path
    test_path = os.path.join(os.path.dirname(__file__), test_dir)
    
    if os.path.exists(test_path):
        suite = loader.discover(test_path, pattern=pattern)
    else:
        print(f"Warning: Test directory {test_path} not found")
        suite = unittest.TestSuite()
    
    return suite


def run_tests(verbosity: int = 2,
             test_pattern: Optional[str] = None,
             fail_fast: bool = False) -> bool:
    """
    Run the test suite.
    
    Args:
        verbosity: Test output verbosity (0-2)
        test_pattern: Pattern to match test files
        fail_fast: Stop on first failure
    
    Returns:
        True if all tests passed, False otherwise
    """
    print("=" * 70)
    print("Gene Curation Pipeline - Test Suite")
    print("=" * 70)
    
    # Discover tests
    pattern = test_pattern or "test_*.py"
    suite = discover_tests(pattern=pattern)
    
    # Count tests
    test_count = suite.countTestCases()
    print(f"Discovered {test_count} tests")
    
    if test_count == 0:
        print("No tests found!")
        return False
    
    print("-" * 70)
    
    # Run tests
    runner = unittest.TextTestRunner(
        verbosity=verbosity,
        failfast=fail_fast,
        buffer=True
    )
    
    start_time = time.time()
    result = runner.run(suite)
    end_time = time.time()
    
    # Print summary
    print("-" * 70)
    print(f"Tests completed in {end_time - start_time:.2f} seconds")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped)}")
    
    if result.failures:
        print("\nFAILURES:")
        for test, traceback in result.failures:
            print(f"  {test}: {traceback}")
    
    if result.errors:
        print("\nERRORS:")
        for test, traceback in result.errors:
            print(f"  {test}: {traceback}")
    
    success = len(result.failures) == 0 and len(result.errors) == 0
    
    if success:
        print("\n✅ All tests passed!")
    else:
        print(f"\n❌ {len(result.failures) + len(result.errors)} test(s) failed")
    
    return success


def run_coverage_tests() -> bool:
    """Run tests with coverage analysis if coverage.py is available."""
    try:
        import coverage
    except ImportError:
        print("Coverage analysis not available (pip install coverage)")
        return run_tests()
    
    print("Running tests with coverage analysis...")
    
    # Start coverage
    cov = coverage.Coverage()
    cov.start()
    
    try:
        # Run tests
        success = run_tests(verbosity=1)
        
        # Stop coverage and generate report
        cov.stop()
        cov.save()
        
        print("\n" + "=" * 70)
        print("COVERAGE REPORT")
        print("=" * 70)
        
        # Generate console report
        cov.report(show_missing=True)
        
        # Generate HTML report if possible
        try:
            cov.html_report(directory="htmlcov")
            print("\nHTML coverage report generated in htmlcov/")
        except:
            pass
        
        return success
        
    except Exception as e:
        print(f"Error running coverage analysis: {e}")
        return False


def run_specific_tests(test_names: List[str]) -> bool:
    """Run specific test classes or methods."""
    suite = unittest.TestSuite()
    
    for test_name in test_names:
        try:
            # Try to load as a test module
            if "." not in test_name:
                test_name = f"gene_curation_pipeline.tests.{test_name}"
            
            parts = test_name.split(".")
            module_name = ".".join(parts[:-1])
            class_name = parts[-1]
            
            module = __import__(module_name, fromlist=[class_name])
            test_class = getattr(module, class_name)
            
            suite.addTest(unittest.TestLoader().loadTestsFromTestCase(test_class))
            
        except (ImportError, AttributeError) as e:
            print(f"Could not load test {test_name}: {e}")
            return False
    
    # Run the specific tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return len(result.failures) == 0 and len(result.errors) == 0


def validate_environment() -> bool:
    """Validate that the test environment is set up correctly."""
    print("Validating test environment...")
    
    # Check required modules
    required_modules = ["unittest", "tempfile", "json", "os", "sys"]
    for module in required_modules:
        try:
            __import__(module)
        except ImportError:
            print(f"❌ Required module {module} not available")
            return False
    
    # Check test files exist
    test_files = [
        "gene_curation_pipeline/tests/test_data_structures.py",
        "gene_curation_pipeline/tests/test_config.py",
        "gene_curation_pipeline/tests/test_codon_integration.py"
    ]
    
    for test_file in test_files:
        if not os.path.exists(test_file):
            print(f"❌ Test file {test_file} not found")
            return False
    
    print("✅ Test environment validation passed")
    return True


def main():
    """Main test runner entry point."""
    parser = argparse.ArgumentParser(description="Run gene curation pipeline tests")
    parser.add_argument("-v", "--verbosity", type=int, choices=[0, 1, 2], default=2,
                       help="Test output verbosity")
    parser.add_argument("-f", "--fail-fast", action="store_true",
                       help="Stop on first failure")
    parser.add_argument("-c", "--coverage", action="store_true",
                       help="Run with coverage analysis")
    parser.add_argument("-p", "--pattern", type=str,
                       help="Pattern to match test files")
    parser.add_argument("-t", "--tests", nargs="+",
                       help="Specific test classes to run")
    parser.add_argument("--validate", action="store_true",
                       help="Only validate test environment")
    
    args = parser.parse_args()
    
    # Validate environment
    if not validate_environment():
        sys.exit(1)
    
    if args.validate:
        print("Test environment is ready!")
        sys.exit(0)
    
    # Run tests
    success = False
    
    if args.tests:
        success = run_specific_tests(args.tests)
    elif args.coverage:
        success = run_coverage_tests()
    else:
        success = run_tests(
            verbosity=args.verbosity,
            test_pattern=args.pattern,
            fail_fast=args.fail_fast
        )
    
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()