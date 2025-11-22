# Fit Empirical Data Using Erlang + Exponential Mixture Model (Version 2)

Fits empirical data to a mixture of Erlang and Exponential distributions
using Maximum Likelihood Estimation (MLE), with options for
goodness-of-fit testing via bootstrapping and adaptive search over shape
parameter `K`.

## Usage

``` r
ErlangExp_Fit_v2(empiricaldata, K, ...)
```

## Arguments

- empiricaldata:

  Numeric vector of observed data.

- K:

  Integer initial guess for the Erlang shape parameter.

- ...:

  Optional named arguments to override defaults:

  Alpha

  :   Significance level for goodness-of-fit test. Default = 0.05.

  pvaloption

  :   Goodness-of-fit metric: `"KS"`, `"CvM"`, `"AD"`, or `"NIL"`.
      Default = `"KS"`.

  InitialguessErLam

  :   Initial guess for Erlang rate parameter. Default = 3.

  SmallestK

  :   Logical. If TRUE, searches for smallest acceptable K under alpha.
      Default = FALSE.

  SmallestKValue

  :   Initial smallest K value for search. Default = -1 (auto).

  ShowFigures

  :   Logical. If TRUE, generates diagnostic plots. Default = TRUE.

  FixedK

  :   Logical. If TRUE, fits only the specified K without adaptive
      search. Default = FALSE.

  KWindowSize

  :   Window size around K for adaptive search. Default = 1.

  NumBootstraps

  :   Number of bootstrap samples for p-value estimation. Default =
      round(10 \* 10 / Alpha).

## Value

A list containing:

- Best:

  A list with the best fit parameters, p-value, and goodness-of-fit
  metrics.

- Smallest:

  If `SmallestK = TRUE`, results for the smallest acceptable K.

- AllFits:

  All fit results across the searched K window (unless `FixedK = TRUE`).

## Details

This function fits an Erlang + Exponential mixture model to the input
data, optionally performs bootstrap-based goodness-of-fit tests, and can
generate diagnostic plots of the fit. It supports an adaptive search
over the shape parameter `K` or fitting with fixed `K`.

## See also

[`ErlangExp_Fit_v2_FixedK()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2_FixedK.md),
[`ErlangExp_Fit_v2_Pvalue()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2_Pvalue.md),
[`ErlangExp_NegLogLikelihood()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_NegLogLikelihood.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- rexp(200, rate = 0.5)
fit <- ErlangExp_Fit_v2(data, K = 2, Alpha = 0.05, ShowFigures = TRUE)
print(fit$Best)
} # }
```
