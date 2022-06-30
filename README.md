
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Implementing ranked sparsity methods with sparseR

[![codecov](https://codecov.io/gh/petersonR/sparseR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/petersonR/sparseR)
[![R-CMD-check](https://github.com/petersonR/sparseR/workflows/R-CMD-check/badge.svg)](https://github.com/petersonR/sparseR/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/sparseR)](https://CRAN.R-project.org/package=sparseR)

## Overview

The `sparseR` package implements the following ranked sparsity methods:

1)  Sparsity-ranked lasso
2)  Sparsity-ranked MCP
3)  Sparsity-ranked SCAD
4)  Ranked Bayesian Information Criterion (and corresponding stepwise
    approaches)

Additionally, `sparseR` has many features designed to streamline dealing
with interaction terms and polynomials, including functions for variable
pre-processing, variable selection, post-selection inference, and
post-fit model visualization under ranked sparsity.

The package is currently in beta phase (version 0.0.1), with plans for a
stable release to CRAN in February 2022. A publication detailing ranked
sparsity principles is currently under review, and available upon
request.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("petersonR/sparseR")
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

For more examples and a closer look, check out the [package
website](https://petersonr.github.io/sparseR/).
