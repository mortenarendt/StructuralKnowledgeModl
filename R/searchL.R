#' Find parameter for Soft thressholding
#' 
#' This functions searches penalty to match sum(abs(u)) contraint via a griddy split half algorithm
#'@param u a vector or matrix and 
#'@param sumabs the required contraint
#'@param niter=100 (default) number of iterations

searchL <- function(u, sumabs,niter = 100) {
  if (norm(as.matrix(u),'E') == 0 || sum(abs(u/norm(as.matrix(u),'E'))) <= sumabs) 
    return(0)
  lam1 <- 0
  lam2 <- max(abs(u))
  iter <- 1
  while (iter < niter) {
    su <- softth(u, (lam1 + lam2)/2)
    if (sum(abs(su/norm(as.matrix(su),'E'))) < sumabs) {
      lam2 <- (lam1 + lam2)/2
    }
    else {
      lam1 <- (lam1 + lam2)/2
    }
    if ((lam2 - lam1) < 1e-06) 
      return((lam1 + lam2)/2)
    iter <- iter + 1
  }
  return((lam1 + lam2)/2)
}
