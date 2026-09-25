# Changelog

## GenErlangFit 0.0.1

*Released 2026-09-24*

#### Bug fixes

- `Erlang_Fit_v2_Pvalue` (`ErlangFuncs.R`): corrected computation of the
  KS and CvM test statistics in sub-sections 1–4 of the function. Prior
  computations did not match base R (`ks.test`) or standard closed-form
  statistical formulations. Statistics are now computed following the
  standard formulations described in D’Agostino & Stephens (1986),
  *Goodness-of-Fit Techniques*.
- `ErlangExp_Fit_v2_Pvalue` (`ErlangExpFuncs.R`): corrected computation
  of the KS and CvM test statistics in sub-sections 1–4 of the function,
  mirroring the fix above, and verified against the same reference.

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
- Declared `ggplot2` and `goftest` as formal package dependencies
  (`Suggests`) in `DESCRIPTION`, reflecting their use in the test suite.

#### Known issues (documented, to be addressed in future patch release)

- `fit$Best$Loglikelihood` (Erlang) vs `fit$Best$LogLikelihood`
  (Erlang-Exp): inconsistent field capitalization across model branches.
  Does not affect core fit or goodness-of-fit functionality.
- Diagnostic plot titles use the literal Unicode character λ, causing a
  locale-dependent rendering failure (`mbcsToSbcs` conversion error)
  under `R CMD check` / non-UTF-8 locales. Does not affect core fit or
  goodness-of-fit functionality; affects only diagnostic plot rendering
  on certain systems.
- Internal diagnostic plotting code uses the deprecated `ggplot2`
  `..density..` syntax. Functional but produces deprecation warnings;
  does not affect core fit or goodness-of-fit functionality.
