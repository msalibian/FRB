Fast and Robust Bootstrap
================
Matias Salibian
2022-05-12

## Fast and Robust Bootstrap

This package implements the Fast and Robust Bootstrap as proposed in
[Salibian-Barrera and Zamar
(2002)](http://dx.doi.org/10.1214/aos/1021379865), and
[Salibian-Barrera, M., Van Aels, S. and Willems, G.
(2008)](http://dx.doi.org/10.1007/s10260-007-0048-6) for robust
regression estimators (MM-estimators) computed with `robustbase::lmrob`.

To install it use the following commands (assuming that you have the
`devtools` package from [CRAN](https://cran.r-project.org) already
installed):

``` r
devtools::install_github("msalibian/FRB")
```

To use it (after installation), simply call `frb` on an `lmrob` object
as computed by `robustbase::lmrob`. Here’s an example:

``` r
library(robustbase)
library(FRB)
a <- lmrob(LNOx ~ LNOxEm + sqrtWS, data=NOxEmissions)
set.seed(123)
tmp <- frb(lmrob.object=a, nboot=1000, return.coef=FALSE)
```

If the argument `return.coef` is set to `FALSE`, then `frb` returns the
estimated covariance matrix of the robust regression estimators. For
example, the estimated standard errors for each parameter estimate are

``` r
sqrt(diag(tmp$var))
```

    ## [1] 0.054340731 0.007633753 0.013364467

We can compare them with the estimated standard errors given by the
usual asyptotic approximation:

``` r
sqrt(diag(summary(a)$cov))
```

    ## (Intercept)      LNOxEm      sqrtWS 
    ## 0.054256788 0.007482346 0.013222502

If the argument `return.coef` is `TRUE` then the returned list has an
element `$coef` that contains a matrix with `nboot` rows with a
bootstrapped regression coefficient in each of them.

The argument `return.indices` indicates whether to return the matrix of
the indices of the bootstrap samples.

The argument `centered` indicates whether the returned bootstrap
regression coefficients should be centered
(*β̂*<sup>\*</sup> − *β̂*<sub>*n*</sub>) or not. For example:

``` r
# by default the bootstrapped estimators are centered
set.seed(123)
tmp <- frb(lmrob.object=a, nboot=1000, return.coef=TRUE)

# we now add the regression estimator to each bootstrapped one
tmp2 <- scale(tmp$coef, scale=FALSE, center=-coef(a))

# use the argument centered=FALSE to return the non-centered 
# bootstrapped estimators instead
set.seed(123)
tmp3 <- frb(lmrob.object=a, nboot=1000, return.coef=TRUE, centered=FALSE)

# and check that the results are the same
all.equal(tmp2, tmp3$coef)
```

    ## [1] TRUE
