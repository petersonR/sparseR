
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Implementing ranked sparsity methods with sparseR

[![codecov](https://codecov.io/gh/petersonR/sparseR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/petersonR/sparseR)
[![R-CMD-check](https://github.com/petersonR/sparseR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/petersonR/sparseR/actions/workflows/R-CMD-check.yaml)
[![CRAN
status](https://www.r-pkg.org/badges/version/sparseR)](https://CRAN.R-project.org/package=sparseR)

## What is ranked sparsity?

The ranked sparsity methods such as the sparsity-ranked lasso (SRL) have
been developed for model selection and estimation in the presence of
interactions and polynomials (Peterson & Cavanaugh
2022)\[<https://doi.org/10.1007/s10182-021-00431-7>\]. The main idea is
that an algorithm should be more skeptical of higher-order polynomials
and interactions a priori compared to main effects, by a predetermined
amount.

## Package overview

The `sparseR` package implements ranked-sparsity-based versions of the
lasso, elastic net, MCP, and SCAD. We also provide a (preliminary)
version of an sparsity-ranked extension to Bayesian Information
Criterion (and corresponding stepwise approaches)

Additionally, `sparseR` has many features designed to streamline dealing
with interaction terms and polynomials, including functions for variable
pre-processing, variable selection, post-selection inference, and
post-fit model visualization under ranked sparsity.

## Installation

``` r

## Via GitHub: 
# install.packages("devtools")
devtools::install_github("petersonR/sparseR")

# or via CRAN
install.packages("sparseR")
```

## Example

``` r
library(sparseR)
```

``` r
data(iris)
set.seed(1321)

srl <- sparseR(Sepal.Width ~ ., data = iris, k = 1, seed = 1)
srl
#> 
#> Model summary @ min CV:
#> -----------------------------------------------------
#>   lasso-penalized linear regression with n=150, p=18
#>   (At lambda=0.0015):
#>     Nonzero coefficients: 10
#>     Cross-validation error (deviance): 0.07
#>     R-squared: 0.62
#>     Signal-to-noise ratio: 1.64
#>     Scale estimate (sigma): 0.267
#> 
#>   SR information:
#>              Vartype Total Selected Saturation Penalty
#>          Main effect     6        4      0.667    2.45
#>  Order 1 interaction    12        6      0.500    3.46
#> 
#> 
#> Model summary @ CV1se:
#> -----------------------------------------------------
#>   lasso-penalized linear regression with n=150, p=18
#>   (At lambda=0.0070):
#>     Nonzero coefficients: 7
#>     Cross-validation error (deviance): 0.08
#>     R-squared: 0.57
#>     Signal-to-noise ratio: 1.33
#>     Scale estimate (sigma): 0.285
#> 
#>   SR information:
#>              Vartype Total Selected Saturation Penalty
#>          Main effect     6        3      0.500    2.45
#>  Order 1 interaction    12        4      0.333    3.46
```

For more examples and a closer look at how to use this package, check
out the [package website](https://petersonr.github.io/sparseR/).

Many thanks to the authors and maintainers of
[`ncvreg`](https://github.com/pbreheny/ncvreg) and
[`recipes`](https://recipes.tidymodels.org/).
