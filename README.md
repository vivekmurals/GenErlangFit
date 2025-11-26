
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GenErlangFit <a href="https://vivekmurals.github.io/GenErlangFit"><img src="man/figures/logo.svg" align="right" height="139" alt="GenErlangFit Logo"></a>

<!-- badges: start -->

<!-- badges: end -->

GenErlangFit is an R package for fitting Erlang and Erlang-Exponential
models.

## Installation

You can install the development version of GenErlangFit from
[GitHub](https://github.com/vivekmurals/GenErlangFit) with:

``` r
install.packages("devtools") # If not already installed.
devtools::install_github("vivekmurals/GenErlangFit", build_vignettes = TRUE)
```

## Documentation

The package documentation is available online at
<https://vivekmurals.github.io/GenErlangFit/>.

The package documentation is also available after installation. For a
quick introduction, see the Get started page. This can also be accessed
from an R session with:

    vignette("GettingStarted", package = "GenErlangFit")

## Data Format and Quick Start

The primary function `GenErlang_Fit` requires at minimum a single
numeric vector of observed time-to-event values where all entries must
be positive. In this documentation, we refer to this vector as
`empiricaldata`.

Once your data is ready, fitting an Erlang or Erlang-Exponential model
is straightforward:

``` r
library(GenErlangFit)
library(ggplot2)

# Default: Automatically fits both Erlang and Erlang-Exponential models,
# performs goodness-of-fit, and outputs the simplest statistically-acceptable model
result <- GenErlang_Fit(empiricaldata)
print(result)

# Fit only the Erlang model and performs goodness-of-fit
result_erlang <- GenErlang_Fit('Erlang', empiricaldata)
print(result_erlang)

# Fit only the Erlang-Exponential model and performs goodness-of-fit 
# (This option requires specifying an initial guess for K)
result_erlangexp <- GenErlang_Fit('ErlangExp', empiricaldata, K = 3)
print(result_erlangexp)
```
