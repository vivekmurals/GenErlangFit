# Bootstrap-Based Goodness-of-Fit for Erlang Distribution

Computes a bootstrap-based p-value for the Erlang model fit using
Kolmogorov–Smirnov, Cramér–von Mises, or Anderson–Darling statistics.

## Usage

``` r
Erlang_Fit_v2_Pvalue(
  empiricaldata,
  k_star,
  lambda_star,
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

  Numeric, fitted Erlang shape parameter.

- lambda_star:

  Numeric, fitted Erlang scale parameter.

- s:

  Integer, number of observations. Defaults to `length(empiricaldata)`.

- n:

  Integer, number of bootstrap samples. Default = 1000.

- alpha:

  Numeric, significance level. Default = 0.05.

- pvaloption:

  Character. Choice of test statistic: `"KS"`, `"CvM"`, or `"AD"`.

- ShowFigures:

  Logical. If TRUE, produces optional plots (currently suppressed).

## Value

A list containing:

- p_value:

  Bootstrap-estimated p-value.

- q_value:

  Binary indicator (1 = fail to reject, 0 = reject).

- metric_star:

  Observed test statistic.

- sample_stats:

  Bootstrap sample statistics.

## Details

The function compares the empirical CDF of the observed data to the
theoretical Erlang CDF and computes a test statistic. Bootstrap samples
from the fitted Erlang model are used to estimate the null distribution
and compute the corresponding p-value.

## See also

[`Erlang_Fit_v2()`](https://vivekmurals.github.io/GenErlangFit/reference/Erlang_Fit_v2.md),
[`Erlang_NegLogLikelihood()`](https://vivekmurals.github.io/GenErlangFit/reference/Erlang_NegLogLikelihood.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- rexp(200, rate = 0.5)
fit <- Erlang_Fit_v2(data, ShowFigures = FALSE)
pvals <- Erlang_Fit_v2_Pvalue(data, fit$Best$K_star, fit$Best$Lambda_star)
} # }
```
