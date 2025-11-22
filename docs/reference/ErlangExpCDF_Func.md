# Compute CDF of Erlang + Exponential Mixture Distribution

Calculates the cumulative distribution function (CDF) values for the
Erlang plus Exponential mixture distribution at specified data points
using numerical integration.

## Usage

``` r
ErlangExpCDF_Func(params, DatX, interval = 0.01)
```

## Arguments

- params:

  Numeric vector of length 3 containing the parameters:

  - Erlang shape parameter `K` (integer)

  - Erlang rate parameter `lambda`

  - Exponential rate parameter `lambda`

- DatX:

  Numeric vector of data points where CDF is evaluated.

- interval:

  Numeric step size for numerical integration grid (default 0.01).

## Value

Numeric vector of CDF values at `DatX`.

## Details

The function computes the PDF values on a grid using `ErlangExp_Func`
and numerically integrates them using the trapezoidal rule to obtain the
CDF. The results are then interpolated at the requested points.

## See also

[`ErlangExp_Func()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Func.md),
[`ErlangExp_Fit_v2_Pvalue()`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2_Pvalue.md)

## Examples

``` r
if (FALSE) { # \dontrun{
params <- c(3, 2, 1)
x <- seq(0, 10, by = 0.1)
cdf_vals <- ErlangExpCDF_Func(params, x)
plot(x, cdf_vals, type = "l")
} # }
```
