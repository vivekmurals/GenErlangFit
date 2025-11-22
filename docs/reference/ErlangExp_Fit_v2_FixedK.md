# Fit Erlang + Exponential Mixture Model with Fixed Shape Parameter K

Fits empirical data to an Erlang + Exponential mixture model using
Maximum Likelihood Estimation (MLE), with a fixed Erlang shape parameter
`K`. Optionally computes bootstrap-based goodness-of-fit p-values.

## Usage

``` r
ErlangExp_Fit_v2_FixedK(empiricaldata, FixedKValue, ...)
```

## Arguments

- empiricaldata:

  Numeric vector of observed data.

- FixedKValue:

  Integer specifying the fixed Erlang shape parameter `K`.

- ...:

  Optional named arguments to override defaults:

  Alpha

  :   Significance level for goodness-of-fit test. Default = 0.05.

  pvaloption

  :   Goodness-of-fit metric: `"KS"`, `"CvM"`, `"AD"`, or `"NIL"` (skip
      p-value). Default = `"KS"`.

  InitialguessErlam

  :   Initial guess for Erlang rate parameter. Default = 3.

  SmallestK

  :   Logical. Not used in this function but included for consistency.
      Default = FALSE.

  ShowFigures

  :   Logical. If TRUE, generates diagnostic plots. Default = TRUE.

  NumBootstraps

  :   Number of bootstrap samples for p-value estimation. Default =
      round(10 \* 10 / Alpha).

## Value

A list containing:

- K_star:

  Estimated Erlang shape parameter (fixed input).

- ErlangLambda_star:

  Estimated Erlang rate parameter.

- ExpLambda_star:

  Estimated Exponential rate parameter.

- P_star:

  Bootstrap p-value for goodness-of-fit (if computed).

- Q_Value:

  Bootstrap q-value (if computed).

- LogLikelihood:

  Log-likelihood of the fitted model.

- metric_star:

  Goodness-of-fit test statistic (if computed).

- samplestats_star:

  Bootstrap sample statistics (if computed).

## Details

This function performs MLE to estimate the Erlang and Exponential rate
parameters, fixing the Erlang shape parameter `K`. Goodness-of-fit
testing via bootstrap is optional and controlled by the `pvaloption`
argument.

## See also

[`ErlangExp_Fit_v2()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2.md),
[`ErlangExp_Fit_v2_Pvalue()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2_Pvalue.md),
[`ErlangExp_NegLogLikelihood()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_NegLogLikelihood.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- rexp(200, rate = 0.5)
fit_fixed <- ErlangExp_Fit_v2_FixedK(data, FixedKValue = 2, Alpha = 0.05)
print(fit_fixed$ErlangLambda_star)
} # }
```
