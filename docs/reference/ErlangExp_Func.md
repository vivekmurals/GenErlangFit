# Erlang + Exponential Mixture PDF and Log-Likelihood

Computes the probability density function (PDF) values of the Erlang
plus Exponential mixture distribution for the given data, and calculates
the total log-likelihood.

## Usage

``` r
ErlangExp_Func(data, ErK, Erlam, Explam)
```

## Arguments

- data:

  Numeric vector of observed data points.

- ErK:

  Integer shape parameter for the Erlang distribution.

- Erlam:

  Numeric rate parameter for the Erlang distribution.

- Explam:

  Numeric rate parameter for the Exponential distribution.

## Value

A list with components:

- Probability:

  Numeric vector of PDF values for each data point.

- Likelihood:

  Sum of log of PDF values (total log-likelihood).

## Details

For each data point, numerically integrates the Erlang + Exponential PDF
using the given parameters. Returns a list containing the vector of PDF
values for each data point and the sum of their log-likelihoods.

## See also

[`ErlangExp_NegLogLikelihood()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_NegLogLikelihood.md),
[`ErlangExp_Fit_v2_FixedK()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2_FixedK.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- rexp(100, rate = 0.5)
res <- ErlangExp_Func(data, ErK = 3, Erlam = 2, Explam = 1)
res$Likelihood
} # }
```
