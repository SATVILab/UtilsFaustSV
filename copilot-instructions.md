# UtilsFaustSV

This is an R package intended for eventual submission to BioConductor.
The package provides additional functions to extract and plot FAUST (Functional Annotation of Unsupervised T-cell clustering) output, allowing users to work with FAUST-identified populations using only a subset of the markers used by FAUST.

## Code standards

### Required before each commit

- Run `devtools::test()` to ensure all tests pass
- Run `devtools::document()` to update documentation
- Run `styler::style_pkg()` to ensure consistent code formatting
- Run `lintr::lint_package()` to check for linting violations

### Development flow

- Build and reload: Use `devtools::load_all()` to reload the package in your R session
- Test: Run `devtools::test()` to execute all unit tests
- Documentation: Use `devtools::document()` to update package documentation
- Coverage: Run `covr::report()` to check test coverage
- Documentation: Run `devtools::document()` to update documentation
- Full check: Run `devtools::check()` to perform a comprehensive package check

## Repository structure

- `R/`: Core R source code
  - `count_get.R`: Functions for extracting counts of FAUST-identified populations
  - `count_plot.R`: Functions for plotting counts of FAUST-identified populations
  - `fcs_write.R`: Functions for saving FAUST-identified populations as FCS files
  - `format.R`: Utility functions for formatting FAUST annotations and marker levels
  - `marker_get.R`: Functions for retrieving marker usage and levels from FAUST output
  - `misc.R`: Miscellaneous utility functions
- `.devcontainer/`: Development container configuration for consistent development environment
- `.github/`: GitHub associated files including workflows and issue templates
- `man/`: Automatically generated documentation files
- `renv/`: R package environment management files
- `scripts/`: Installation and setup scripts for development environment
- `tests/`: Unit tests for the package
  - `testthat/`: Test files using the testthat framework
    - Test data directories: `faustData/`, `fcsTest/`, `gs_test/`
- `_dependencies.R`: Explicitly listed dependencies for `renv` to pick up via `implicit` dependencies
- `.gitignore`: Git ignore file to exclude unnecessary files from version control
- `.Rbuildignore`: R build ignore file to exclude unnecessary files from package builds
- `DESCRIPTION`: Package metadata file
- `LICENSE`: License file for the package
- `LICENSE.md`: License file in Markdown format
- `README.md`: Package overview and usage documentation
- `README.Rmd`: R Markdown source for README.md
- `renv.lock`: Lock file for the `renv` package management system, capturing the state of the R package environment

## Key Guidelines

1. Begin each internal function with a `.`
2. Use `.debug` as the parameter name for debugging flags in internal functions
3. Use `.debug_msg()` for debug messages, which takes a boolean, a message and an optional value
4. Add unit tests using `testthat` for all new functionality
5. Validate inputs and provide meaningful error messages
6. Explicitly refer to all packages used, rather than using `@import` or `@importFrom`, with the exception of `ggplot2` functions, the `flowCore::exprs` function.
7. Use `@export` for functions that should be available to users
8. When running code from the project, you must always have as your working directory the root of the project, i.e. the directory containing the `DESCRIPTION` file. This is especially important when the project uses `renv`, as otherwise the `.Rprofile` will not be sourced and the package environment will not be set up correctly.
8. Never update `.Rd` files manually; use `devtools::document()` to regenerate them.
9. Documentation in `.R` files for parameters must always have the format `@param param_name <type_of_input> <info>`.
  - For example, `@param project_path character Path to FAUST project directory.`,
  - Of course, the information can be longer, and where there are multiple options, it can be mentioned:
    - For example, `logical or character`, or `"always", "never" or "automatic"`.
  - If there is a default specified, it should be stated at the end in the documentation of that parameter, e.g. `Default is "automatic".`.
  - Where longer parameter descriptions are required, usually you should describe overall what the parameter is for, and then list the options with a short description of each.
  - If there are really nitty-gritty details, then you can use `@details` to provide more information and refer to it in the parameter description.
10. Always make sure there is no unnecessary trailing whitespace, including for blank lines.
11. Always make sure that any edited files have a final newline at the end of the file.
12. Tests written should not simply test whether it errors out correctly, but should also test that the output is as expected.
13. Never use `return` for the last line of a function, but only when you want to return early from a function.
14. For testing, never add source commands at the top to source files in the `R` folder, as these are automatically sourced (effectively) by the `testthat` package.
15. The tests need to pass on Mac, Windows and Ubuntu, so be aware of that (e.g. file path availability and specification (forward or backward slashes, root directories, etc.), case sensitivity, etc.). You do not run the tests on all three platforms, but you should think about them passing on all three platforms.

## FAUST-specific guidelines

16. When working with FAUST population definitions, use the tilde notation (e.g., "CD3~2~2~") for internal processing and convert to/from plus/minus notation (e.g., "CD3+") for user-facing interfaces.
17. Population definitions can be specified as:
    - Named character vectors (e.g., `c("CD3" = 2, "CD4" = 2)`)
    - Lists of named character vectors for multiple populations
    - Single character strings with full FAUST annotations
18. Always validate that the specified markers exist in the FAUST output before processing.
19. When working with GatingSets, ensure proper handling of sample names and indices.
20. FCS file outputs should maintain the original marker names and expression values unless explicitly transformed.