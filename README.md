# Fast and Robust Bootstrap

This package implements the Fast and Robust Bootstrap as proposed in 
[Salibian-Barrera and Zamar (2002)](http://dx.doi.org/10.1214/aos/1021379865), and
[Salibian-Barrera, M., Van Aels, S. and Willems, G. (2008)](http://dx.doi.org/10.1007/s10260-007-0048-6) for robust regression estimators (MM-estimators) computed with 
`robustbase::lmrob`. 

To install it use the following commands (assuming that you have
the `devtools` package from [CRAN](https://cran.r-project.org) 
already installed): 
```
library(devtools)
install_github("msalibian/FRB")
```
To use it (after installation), simply call `frb` on an `lmrob` object as computed 
by `robustbase::lmrob`. Here's an example:
```
library(robustbase)
library(FRB)
a <- lmrob(LNOx ~ LNOxEm + sqrtWS, data=NOxEmissions)
tmp <- frb(lmrob.object=a, nboot=1000, return.coef=FALSE)
# Estimated SE's for estimated regression coefficients  
sqrt(diag(tmp))
# [1] 0.056422169 0.007782671 0.012662991
# compare with SE's based on the asymptotic approximation
sqrt(diag(summary(a)$cov))
# (Intercept)      LNOxEm      sqrtWS 
# 0.054256788 0.007482346 0.013222502 
```
