# Negative Log-Likelihood for Erlang + Exponential Mixture

Computes the negative log-likelihood of the Erlang + Exponential mixture
distribution given the Erlang rate parameter, data, and Erlang shape
parameter.

## Usage

``` r
ErlangExp_NegLogLikelihood(ErlangLambda, data, K)
```

## Arguments

- ErlangLambda:

  Numeric rate parameter for the Erlang component.

- data:

  Numeric vector of observed data.

- K:

  Integer shape parameter for the Erlang distribution.

## Value

Numeric value of the negative log-likelihood.

## Details

Calculates the exponential rate parameter from the mean constraint and
evaluates the negative log-likelihood using the Erlang + Exponential PDF
implemented in `ErlangExp_Func`. Returns `Inf` if the exponential
parameter is invalid (e.g., negative or infinite) to penalize invalid
parameter sets.

## See also

[`ErlangExp_Func()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Func.md),
[`ErlangExp_Fit_v2_FixedK()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2_FixedK.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- rexp(100, rate = 0.5)
negLL <- ErlangExp_NegLogLikelihood(3, data, 2)
} # }
```
