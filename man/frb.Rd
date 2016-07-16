\name{frb}
\alias{frb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to run the Fast and Robust Bootstrap  }
\description{
This function runs a Fast and Robust Bootstrap for
robust regression estimators (MM) as computed by
\code{robustbase::lmrob}. 
}
\usage{
frb(lmrob.object, nboot = 1000, return.coef = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lmrob.object}{ An object of class \code{lmrob} as
  returned by \code{robustbase::lmrob} }
  \item{nboot}{ An integer number. The number of bootstrap samples to be used.}
  \item{return.coef}{ A logical value. If \code{TRUE} then the function returns
  a matrix of 
  \code{nboot} by \code{p}, with the i-th row corresponding to 
  the i-th  bootstrapped 
  regression estimators. If \code{FALSE} (default) then the function
  returns the estimated covariance matrix (\code{p} by \code{p}) for the
  MM-regression estimator.}
}
\details{
The fast and robust bootstrap as described in 
Salibian-Barrera, M. and Zamar, R.H. (2002), and 
Salibian-Barrera, M., Van Aels, S. and Willems, G. (2008). 
}
\value{ 
Depending on the value of the argument \code{return.coef}, either 
the estimated covariance matrix for the MM-regression estimator, or
a matrix containing the \code{nboot} re-computed robust regression
estimators.
}
\references{
Salibian-Barrera, M. and Zamar, R.H. (2002). Bootstrapping robust estimates of regression. The Annals of Statistics, 30, 556-582.
http://dx.doi.org/10.1214/aos/1021379865

Salibian-Barrera, M., Van Aels, S. and Willems, G. (2008). Fast and robust bootstrap. Statistical Methods and Applications 17, 41-71. 
http://dx.doi.org/10.1007/s10260-007-0048-6
}
\author{Matias Salibian-Barrera, matias@stat.ubc.ca}
\note{See the github repository \link{https://github.com/msalibian/FRB}
}

\keyword{robustness}
\keyword{robust regression}
\keyword{bootstrap}% __ONLY ONE__ keyword per line

