# Changelog

## GenErlangFit 0.0.1

*Released 2026-09-25*

#### Bug fixes

- The Kolmogorov-Smirnov and Cramer-von Mises goodness-of-fit statistics
  (`Erlang_Fit_v2_Pvalue` in `ErlangFuncs.R`; `ErlangExp_Fit_v2_Pvalue`
  in `ErlangExpFuncs.R`) have been revised to correctly implement their
  standard closed-form definitions, following D’Agostino & Stephens
  (1986), *Goodness-of-Fit Techniques*. Independent unit tests
  validating these statistics against base R (`ks.test`) and the
  `goftest` package (`cvm.test`, `ad.test`) have been added to the test
  suite to confirm correctness and prevent regression.
- Diagnostic plot titles previously used the literal Unicode character
  λ, causing a locale-dependent rendering failure (`mbcsToSbcs`
  conversion error) under `R CMD check` on non-UTF-8 locales. Replaced
  with the word “lambda” throughout (`ErlangFuncs.R`,
  `ErlangExpFuncs.R`).
- Removed non-ASCII characters (en dashes, curly apostrophes, accented
  characters, multiplication signs) from source files and roxygen
  documentation (`ErlangFuncs.R`, `ErlangExpFuncs.R`, `MainDriver.R`,
  `DESCRIPTION`), which triggered an `R CMD check` warning.
- Internal diagnostic plotting code used the deprecated `ggplot2`
  `..density..` syntax; updated to `after_stat(density)`.
- Corrected missing package imports: added explicit `@importFrom`
  declarations for functions from `stats` and `utils`, and formally
  declared `ggplot2` as an `Imports` dependency in `DESCRIPTION`
  (previously undeclared, causing an `R CMD check` NOTE for undefined
  global functions).

#### Testing

- Added a formal `testthat` unit test suite
  (`tests/testthat/test-core-model-fitting.R`), comprising 9 tests:
  - 4 tests for the Erlang model: log-likelihood, KS, AD, and CvM
    statistics, each cross-validated against independent base R /
    `goftest` implementations.
  - 1 foundational test validating the Erlang-Exponential CDF helper
    function (`ErlangExpCDF_Func`) against an independent Monte Carlo
    simulation, to avoid circular validation in the tests below.
  - 4 tests for the Erlang-Exponential model: log-likelihood, KS, AD,
    and CvM statistics, cross-validated against the same external
    references.
- Added a second unit test suite
  (`tests/testthat/test-model-selection-logic.R`), comprising 15 tests
  covering model-selection and orchestration logic not exercised by the
  core statistical tests above:
  - 2 tests validating the “Smallest K” goodness-of-fit selection logic
    (Erlang and Erlang-Exponential models), confirming the selected K
    passes goodness-of-fit while the next-smaller K correctly fails.
  - 1 test validating the adaptive window K-search peak-finding logic
    for the Erlang-Exponential model, confirming the selected K is a
    genuine local maximum in log-likelihood.
  - 1 test validating the K-clamping formula used in the “quick fit all
    models” convenience mode.
- Total test suite: 24 tests, all passing with 0 warnings.
- Package test coverage (via `covr`) increased from 53.9% to 73.1%
  overall; core statistical fitting files at 78.4% (`ErlangFuncs.R`) and
  82.9% (`ErlangExpFuncs.R`).
- Declared `ggplot2` and `goftest` as formal package dependencies
  (`Suggests`) in `DESCRIPTION`, reflecting their use in the test suite.

#### Package maintenance

- Corrected a filename typo in `.Rbuildignore` that caused the RStudio
  project file to be incorrectly excluded/included in build checks.
- Removed a stray `.RData` autosave file from the bundled Shiny app
  directory.
- Fixed NEWS.md formatting (unparseable duplicate header) flagged by
  `R CMD check`.
- `R CMD check` now returns 0 errors and 0 warnings, with 1 remaining
  note (unable to verify system clock against CRAN time server), which
  is an environmental check unrelated to package code.

#### Known issues

- None
