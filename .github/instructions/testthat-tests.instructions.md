---
applyTo: "tests/testthat/test-*.R"
---

## testthat Test Requirements

When writing or modifying testthat tests, please follow these guidelines:

1. **Test structure and organization**
   - Use `test_that()` for each test case with descriptive names
   - Group related tests together in the same file
   - File names must follow the pattern: `test-<feature>.R`
   - Never add `source()` commands to load R files - testthat automatically sources them

2. **Test quality and coverage**
   - Tests must validate actual output, not just error handling
   - Use specific assertions that check expected values, not just that code runs
   - Test both success cases and error cases
   - Include edge cases in your tests

3. **Assertions**
   - Use appropriate `expect_*()` functions:
     - `expect_identical()` for exact matches
     - `expect_equal()` for near-equality (numeric)
     - `expect_error()` for expected errors
     - `expect_warning()` for expected warnings
     - `expect_true()` / `expect_false()` for boolean tests
   - Always include descriptive failure messages when helpful

4. **Cross-platform compatibility**
   - Tests must pass on Mac, Windows, and Ubuntu
   - Be aware of:
     - File path differences (forward vs backward slashes)
     - Case sensitivity differences
     - Root directory specifications
     - Line ending differences
   - Use `file.path()` for constructing paths
   - Avoid hardcoded absolute paths

5. **Test data**
   - Use test data from existing directories in `tests/testthat/`:
     - `faustData/`: FAUST output test data
     - `fcsTest/`: FCS file test data
     - `gs_test/`: GatingSet test data
   - These directories exist and contain the necessary test fixtures
   - Keep test data minimal and focused
   - Do not commit large test files without review

6. **Code style in tests**
   - Remove all unnecessary trailing whitespace
   - Ensure files end with a final newline
   - Follow the same style conventions as package code

7. **Working directory**
   - When running tests, working directory must be the project root (directory containing `DESCRIPTION`)
   - This is critical for `renv` to work correctly
   - Tests should not change the working directory
