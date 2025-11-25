
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

## Example

Here is a simple example showing basic usage:

``` r
library(GenErlangFit)
result <- GenErlang_Fit(data)
print(result)
```
