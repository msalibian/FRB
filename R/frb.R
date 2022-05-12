#' Fast and Robust Bootstrap
#'
#' This function computes the Fast and Robust Bootstrap for linear models
#' using etc. etc. It relies on the function \code{lmrob} of 
#' \code{package:robustbase}
#'
#' This function computes the Fast and Robust Bootstrap for linear models
#' using etc. etc. It relies on the function \code{lmrob} of 
#' \code{package:robustbase}
#'
#' @param lmrob.object a robust linear regression fit as returned by \code{lmrob} of \code{package:robustbase}
#' @param nboot an integer indicating the number of bootstrap samples to use
#' @param return.coef a logical (boolean) value indicating whether to return a matrix containing
#' the \code{nboot} vectors of bootstrapped regression coefficient estimates. 
#' @param return.indices a logical (boolean) value indicating whether to return a matrix containing
#' the vectors of indices corresponding to each bootstrap sample used. 
#' @param centered a logical value indicating whether to return the centered bootstrapped
#' regression coefficients (centered around the MM regression estimator). Defaults to \code{TRUE}.
#'
#'
#' @return A list with the following possible components:
#' \item{coef}{If the argument \code{return.coef == TRUE}, a matrix with 
#' \code{nboot} rows containing the bootstrapped robust regression estimates}
#' \item{var}{The sample covariance matrix of the bootstrapped regression estimates. This
#' is an estimate for the asymptotic covariance matrix of the robust regression estimator.}
#' \item{indices}{If the argument \code{return.indices == TRUE}, a matrix with 
#' \code{nboot} rows containing the indices of the bootstrap samples used.}
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}
#' @references Salibian-Barrera M, Zamar RH, Bootstrapping robust estimates of regression,
#' Annals of Statistics, (2002) 30:556-582. \doi{10.1214/aos/1021379865}
#' Salibian-Barrera M, Van Aelst S, Willerns G., Fast and robust bootstrap, 
#' Statistical Methods and Applications, (2008) 17:41-71. \doi{10.1007/s10260-007-0048-6}
#' 
#' @seealso \code{package:robustbase}
#'
#' @examples
#' library(robustbase)
#' a <- lmrob(LNOx ~ LNOxEm + sqrtWS, data=NOxEmissions)
#' set.seed(123)
#' tmp <- frb(lmrob.object=a, nboot=1000, return.coef=FALSE)
#' sqrt(diag(tmp$var))
#' # compare with SE's based on the asymptotic approximation
#' sqrt(diag(summary(a)$cov))
#
#' @export
frb <- function(lmrob.object, nboot=1000, return.coef = FALSE, 
                return.indices = FALSE, centered=TRUE) {
  lmrob.Chi <- Mchi
  lmrob.Psi <- Mpsi
  co <- lmrob.object$control
  tuning.psi <- co$tuning.psi
  tuning.chi <- co$tuning.chi
  beta.s <- lmrob.object$init.S$coefficients
  beta.mm <- coef(lmrob.object)
  re.mm <- as.vector( residuals(lmrob.object) )
  uu <- model.frame(lmrob.object)
  yy <- as.vector( model.extract(uu, 'response') )
  xx <- model.matrix(lmrob.object)
  n <- (dd <- attr(xx, 'dim'))[1]
  p <- dd[2]
  attributes(xx) <- NULL
  xx <- matrix(xx, n, p )
  re.s <- as.vector( yy - xx %*% beta.s )
  scale <- lmrob.object$scale
  # w.p <- Psi'(r/sigma)
  w.p <- lmrob.Psi(x = re.mm / scale, psi='bisquare', cc=tuning.psi, deriv=1) 
  # a more efficient way of doing this?
  a <- t(xx) %*% diag(w.p) %*% xx 
  # w <- Psi(r/sigma)/r
  w <- lmrob.Psi(x= re.mm / scale, psi='bisquare', cc=tuning.psi, deriv=0) / re.mm 
  # a more efficient way of doing this?
  b <- t(xx) %*% diag(w) %*% xx 
  # w.pp <- Psi'(r/sigma)*r/sigma
  w.pp <- w.p * re.mm / scale  
  d <- as.vector( t(xx) %*% w.pp ) 
  # e <- \sum Chi'(r/sigma)*r/sigma
  e <- sum( lmrob.Chi(re.s / scale, psi='bisquare', cc=tuning.chi, deriv=1) * re.s / scale ) 
  d <- d / 2 * (n - p) * scale / e;
  x2 <- solve(a) 
  x3 <- x2 %*% b * scale
  v2 <- x2 %*% d # / scale
  # correction matrix is in x3
  # correction vector is in v2
  ss <- lmrob.Chi( re.s/scale, psi='bisquare', cc=tuning.chi, deriv=0 ) 
  # ss <- Chi(re.s/scale)
  boot.beta <- matrix(0, nboot, p)
  boot.indices <- matrix(0, nboot, n)
  a <- .C('R_frb', as.double(xx), as.double(yy), as.double(w), as.integer(n), 
          as.integer(p), as.double(beta.mm), as.double(scale), 
          as.double(ss), bb=as.double(boot.beta), as.integer(nboot),
          as.double(x3), as.double(v2), 
          bind = as.integer(boot.indices), PACKAGE="FRB")
  ab <- matrix(a$bb, nboot, p)
  if(!centered) ab <- scale(ab, scale=FALSE, center = -beta.mm)
  ib <- matrix(a$bind, nboot, n)
  tmp <- list(coef=ab, var=var(ab), indices = ib)
  if(return.coef) 
    tmp$var <- NULL 
  else tmp$coef <- NULL
  if(!return.indices) tmp$indices <- NULL
  return(tmp)
}
