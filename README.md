Fast and Robust Bootstrap
================
Matias Salibian
2016-07-30

Fast and Robust Bootstrap
-------------------------

This package implements the Fast and Robust Bootstrap as proposed in [Salibian-Barrera and Zamar (2002)](http://dx.doi.org/10.1214/aos/1021379865), and [Salibian-Barrera, M., Van Aels, S. and Willems, G. (2008)](http://dx.doi.org/10.1007/s10260-007-0048-6) for robust regression estimators (MM-estimators) computed with `robustbase::lmrob`.

To install it use the following commands (assuming that you have the `devtools` package from [CRAN](https://cran.r-project.org) already installed):

``` r
library(devtools)
install_github("msalibian/FRB")
```

To use it (after installation), simply call `frb` on an `lmrob` object as computed by `robustbase::lmrob`. Here's an example:

``` r
library(robustbase)
library(FRB)
a <- lmrob(LNOx ~ LNOxEm + sqrtWS, data=NOxEmissions)
set.seed(123)
tmp <- frb(lmrob.object=a, nboot=1000, return.coef=FALSE)
```

If the argument `return.coef` is set to `FALSE`, then `frb` returns the estimated covariance matrix of the robust regression estimators. For example, the estimated standard errors for each parameter estimate are

``` r
sqrt(diag(tmp))
```

    ## [1] 0.054340731 0.007633753 0.013364467

We can compare them with the estimated standard errors given by the usual asyptotic approximation:

``` r
sqrt(diag(summary(a)$cov))
```

    ## (Intercept)      LNOxEm      sqrtWS 
    ## 0.054256788 0.007482346 0.013222502
