# Bootstrap-based Goodness-of-Fit Test for Erlang + Exponential Fit

Performs a bootstrap hypothesis test to assess goodness-of-fit of
empirical data to an Erlang + Exponential mixture model, using the
specified test statistic (Kolmogorov-Smirnov, Cramér-von Mises, or
Anderson-Darling).

## Usage

``` r
ErlangExp_Fit_v2_Pvalue(
  empiricaldata,
  k_star,
  erlambda_star,
  explambda_star,
  s = length(empiricaldata),
  n = 1000,
  alpha = 0.05,
  pvaloption = "KS",
  ShowFigures = TRUE
)
```

## Arguments

- empiricaldata:

  Numeric vector of observed data.

- k_star:

  Integer Erlang shape parameter (fixed).

- erlambda_star:

  Numeric Erlang rate parameter.

- explambda_star:

  Numeric Exponential rate parameter.

- s:

  Integer sample size for bootstrap samples. Default is length of
  empiricaldata.

- n:

  Integer number of bootstrap samples to generate. Default is 1000.

- alpha:

  Significance level for hypothesis testing. Default is 0.05.

- pvaloption:

  Character specifying goodness-of-fit metric to use: `"KS"`
  (Kolmogorov-Smirnov), `"CvM"` (Cramér-von Mises), or `"AD"`
  (Anderson-Darling). Default is `"KS"`.

- ShowFigures:

  Logical indicating whether to generate diagnostic plots. Default is
  TRUE (currently commented out in code).

## Value

A list containing:

- p_value:

  Bootstrap p-value for the goodness-of-fit test.

- q_value:

  Binary decision indicator: 1 = fail to reject null, 0 = reject null.

- metric_star:

  Observed test statistic value from empirical data.

- sample_stats:

  Vector of test statistics from bootstrap samples.

## Details

This function computes the empirical CDF and theoretical
Erlang+Exponential CDF, calculates the specified goodness-of-fit test
statistic, and then uses bootstrap resampling to estimate the p-value
for the test. It also provides a simple hypothesis test decision based
on the significance level `alpha`.

## See also

[`ErlangExp_Fit_v2()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2.md),
[`ErlangExp_Fit_v2_FixedK()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2_FixedK.md),
[`ErlangExpCDF_Func()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExpCDF_Func.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- rexp(200, rate = 0.5)
fit <- ErlangExp_Fit_v2(data)
pval_res <- ErlangExp_Fit_v2_Pvalue(data, fit$Best$K_star, fit$Best$ErlangLambda_star, fit$Best$ExpLambda_star)
print(pval_res$p_value)
} # }
```
