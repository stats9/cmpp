# cmpp 0.0.2

## Improvements and Fixes

### 🛠️ Fixed Examples
- Corrected runtime issues in the following example functions that previously caused ASAN-triggered errors due to mismatched input lengths:
  - `F_cdf_rcpp()`, `F_cdf_rcpp2()`, `F_cdf_rcpp3()`
  - `f_pdf_rcpp()`, `f_pdf_rcpp2()`, `f_pdf_rcpp3()`
- Updated examples to ensure proper object initialization and matching dimensions.

### 🔒 C++ Safety Enhancements
- Added internal length checks (`Z.size() == Beta.size()`) in:
  - `F_cdf()`, `F_cdf2()`, `F_cdf3()`
  - `f_pdf()`, `f_pdf2()`, `f_pdf3()`
- These checks prevent buffer overflows and improve robustness under memory sanitizers.Stop in 1980. 

### 🧪 Function Logic Refinement
- Improved `Cmpp_CIF()` to resolve inconsistencies in event-level handling and ensure proper matrix manipulation.

### 📦 Package Metadata
- Minor edit to the `Title` field in `DESCRIPTION` to better reflect the functionality.
- All checks (`R CMD check --as-cran`) pass with no ERROR or WARNING.
