# Generalized Erlang and Erlang-Exponential Distribution Fitting

Fits empirical data to Erlang or Erlang-Exponential mixture models using
maximum likelihood estimation. This wrapper function calls the
appropriate specialized fitting routines based on the `mode` parameter.

## Usage

``` r
GenErlang_Fit(mode, empiricaldata = NULL, K = NULL, ...)
```

## Arguments

- mode:

  Character string specifying the model to fit. Options are `"Erlang"`
  or `"ErlangExp"`. Alternatively, if `mode` is numeric with no
  additional arguments, it is interpreted as empirical data, and both
  models are fit.

- empiricaldata:

  Numeric vector of empirical data to be fitted. Required when `mode` is
  `"Erlang"` or `"ErlangExp"`.

- K:

  Integer specifying the initial guess for the shape parameter `K`
  (number of Erlang phases) for `"ErlangExp"` mode. Required if
  `mode = "ErlangExp"`.

- ...:

  Additional optional parameters passed down to the underlying fitting
  functions. For example:

  - `SmallestK` Logical. Whether to search for the smallest acceptable
    `K` satisfying goodness-of-fit.

  - Other model-specific options.

## Value

A list containing fitting results:

- For `"Erlang"` mode, returns the result list from
  [`Erlang_Fit_v2`](https://vivekmurals.github.io/GenErlangFit/reference/Erlang_Fit_v2.md).

- For `"ErlangExp"` mode, returns the result list from
  [`ErlangExp_Fit_v2`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2.md).

- If only empirical data is provided (numeric `mode`), returns a list
  with both Erlang and Erlang-Exp fit results.

## See also

[`Erlang_Fit_v2`](https://vivekmurals.github.io/GenErlangFit/reference/Erlang_Fit_v2.md),
[`ErlangExp_Fit_v2`](https://vivekmurals.github.io/GenErlangFit/reference/ErlangExp_Fit_v2.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data <- rexp(100, rate = 0.1)
GenErlang_Fit("Erlang", data, SmallestK = TRUE)
GenErlang_Fit("ErlangExp", data, K = 3)
GenErlang_Fit(data) # fits both models and returns results
} # }
```
