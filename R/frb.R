
frb <- function(lmrob.object, nboot=1000, return.coef = FALSE) # , seed=33) 
{
  lmrob.Psi <- tukeyPsi1
  lmrob.Chi <- tukeyChi
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
  w.p <- lmrob.Psi(x = re.mm / scale, cc=tuning.psi, deriv=1) 
  # a more efficient way of doing this?
  a <- t(xx) %*% diag(w.p) %*% xx 
  # w <- Psi(r/sigma)/r
  w <- lmrob.Psi(x= re.mm / scale, cc=tuning.psi, deriv=0) / re.mm 
  # a more efficient way of doing this?
  b <- t(xx) %*% diag(w) %*% xx 
  # w.pp <- Psi'(r/sigma)*r/sigma
  w.pp <- w.p * re.mm / scale  
  d <- as.vector( t(xx) %*% w.pp ) 
  # e <- \sum Chi'(r/sigma)*r/sigma
  e <- sum( lmrob.Chi(re.s / scale, cc=tuning.chi, deriv=1) * re.s / scale ) 
  d <- d / 2 * (n - p) * scale / e;
  x2 <- solve(a) 
  x3 <- x2 %*% b * scale
  v2 <- x2 %*% d # / scale
  # correction matrix is in x3
  # correction vector is in v2
  # set.seed(seed)
  ss <- lmrob.Chi( re.s/scale, cc=tuning.chi, deriv=0 ) 
  # ss <- Chi(re.s/scale)
  boot.beta <- matrix(0, nboot, p)
  a <- .C('R_frb', as.double(xx), as.double(yy), as.double(w), as.integer(n), 
          as.integer(p), as.double(beta.mm), as.double(scale), 
          as.double(ss), bb=as.double(boot.beta), as.integer(nboot),
          as.double(x3), as.double(v2), PACKAGE="FRB")$bb
  a <- matrix(a, nboot, p)
  if(return.coef) return(a) 
  else return(var(a))
}


