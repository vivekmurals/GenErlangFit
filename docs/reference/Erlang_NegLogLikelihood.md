# Negative Log-Likelihood for Erlang Distribution

Computes the negative log-likelihood (to be minimized) for a given
Erlang (Gamma) shape parameter `k` and data vector. The scale parameter
is internally estimated as \\lambda = mean(data) / k\\.

## Usage

``` r
Erlang_NegLogLikelihood(k, data)
```

## Arguments

- k:

  Numeric, Erlang shape parameter.

- data:

  Numeric vector of observed data.

## Value

Numeric scalar representing the negative log-likelihood.

## See also

[`Erlang_Fit_v2()`](https://vivekmurals.github.io/GenErlangFit/reference/Erlang_Fit_v2.md),
[`Erlang_Fit_v2_Pvalue()`](https://vivekmurals.github.io/GenErlangFit/reference/Erlang_Fit_v2_Pvalue.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- rexp(200, rate = 0.5)
Erlang_NegLogLikelihood(3, data)
} # }
```
