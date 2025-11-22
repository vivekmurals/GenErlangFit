# Fit Empirical Data Using Erlang Distribution (Version 2)

Fits empirical data to an Erlang (Gamma with integer shape) distribution
using Maximum Likelihood Estimation (MLE), and evaluates goodness-of-fit
using bootstrap-based hypothesis testing.

## Usage

``` r
Erlang_Fit_v2(empiricaldata, ...)
```

## Arguments

- empiricaldata:

  Numeric vector of observed data.

- ...:

  Optional named arguments to override defaults:

  Alpha

  :   Significance level for goodness-of-fit test. Default = 0.05.

  pvaloption

  :   Goodness-of-fit metric: `"KS"`, `"CvM"`, or `"AD"`. Default =
      `"KS"`.

  InitialguessK

  :   Initial guess for Erlang shape parameter. Default = 3.

  SmallestK

  :   Logical. If TRUE, searches for smallest acceptable K under alpha.
      Default = FALSE.

  ShowFigures

  :   Logical. If TRUE, generates plots. Default = TRUE.

  NumBootstraps

  :   Number of bootstrap samples for p-value estimation. Default =
      round(10 \* 10 / Alpha).

## Value

A list containing:

- Best:

  A list with MLE estimates, p-value, log-likelihood, and test metrics.

- Smallest:

  If `SmallestK = TRUE`, results for the smallest acceptable K.

## Details

This function estimates the Erlang parameters (shape `K` and scale
`lambda`) by maximizing the log-likelihood function, and tests whether
the data are consistent with the fitted model via bootstrap-based
Kolmogorov–Smirnov (or other) tests. Optionally, it can also search for
the smallest integer `K` that still passes the chosen significance
threshold (`SmallestK = TRUE`).

## See also

[`Erlang_Fit_v2_Pvalue()`](https://vivekmurals.github.io/GenErlangFit/reference/Erlang_Fit_v2_Pvalue.md),
[`Erlang_NegLogLikelihood()`](https://vivekmurals.github.io/GenErlangFit/reference/Erlang_NegLogLikelihood.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- rexp(200, rate = 0.5)
fit <- Erlang_Fit_v2(data, Alpha = 0.05, pvaloption = "KS", ShowFigures = FALSE)
fit$Best$K_star
} # }
```
