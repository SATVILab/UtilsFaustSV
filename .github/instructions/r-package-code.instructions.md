---
applyTo: "R/*.R"
---

## R Package Code Requirements

When writing or modifying R package source code, please follow these guidelines:

1. **Function naming conventions**
   - Begin each internal function with a `.` (e.g., `.internal_helper_function`)
   - Use descriptive, snake_case names for all functions
   - Export functions intended for users with `@export` roxygen tag

2. **Documentation standards**
   - All parameters must follow the format: `@param param_name <type_of_input> <info>`
   - Example: `@param project_path character Path to FAUST project directory.`
   - For parameters with multiple options, list them: `logical or character` or `"always", "never" or "automatic"`
   - If a default exists, state it at the end: `Default is "automatic".`
   - Use `@details` for nitty-gritty details, not parameter descriptions
   - Never manually edit `.Rd` files; use `devtools::document()` to regenerate them

3. **Package dependencies**
   - Explicitly refer to all packages used with `package::function()` notation
   - Exception: Core `ggplot2` functions (like `ggplot()`, `aes()`, `geom_*()`, etc.) and `flowCore::exprs` can be used without prefix
   - Do NOT use `@import` or `@importFrom` (except for ggplot2 and flowCore::exprs)

4. **Debugging support**
   - Use `.debug` as the parameter name for debugging flags in internal functions
   - Use `.debug_msg(.debug, message, value = NULL)` for debug messages (if this utility function is available)

5. **Code quality**
   - Remove all unnecessary trailing whitespace, including on blank lines
   - Ensure all files end with a final newline
   - Never use `return` for the last line of a function; only use it for early returns
   - Validate all inputs and provide meaningful error messages using `stop()`

6. **FAUST-specific requirements**
   - Use tilde notation (e.g., "CD3~2~2~") for internal processing
   - Convert to/from plus/minus notation (e.g., "CD3+") for user-facing interfaces
   - Population definitions can be:
     - Named character vectors: `c("CD3" = 2, "CD4" = 2)`
     - Lists of named character vectors for multiple populations
     - Single character strings with full FAUST annotations
   - Always validate that specified markers exist in FAUST output before processing
   - When working with GatingSets, ensure proper handling of sample names and indices
   - FCS file outputs must maintain original marker names and expression values unless explicitly transformed
