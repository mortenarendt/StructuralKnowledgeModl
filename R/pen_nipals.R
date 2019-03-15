#' NIPALS algorithm for sparse bilinear one component model of a matrix
#'
#' This functions minimises |x-udv'|2 s.t. sum(abs(u))<sumabs[1], sum(abs(v))<sumabs[2], norm(u,2)=norm(v,2) = 1 with respect to u and v.
#'@param x a n by p dataset
#'@param sumabs 1 by 2 vector ofL1 contraint on u and v respectively
#'@param niter = 100 (default) number of iterations
#'
#'@example
#'
#'x <- matrix(rnorm(200),10,20)
#'results <- pen_nipals(x, sumabs = c(2,3))
#'
#'@export
pen_nipals <- function (x,sumabs = sqrt(c(nrow(x),ncol(x))), niter = 100){
  # get LS initials
  r <- svd(x,nu = 1,nv = 1)
  D <- c(0,0)
  c <- 2
  while ((D[c] - D[c-1])>=0 & c<niter){
    c <- c+1
    # update u
    r$u <- smr(x,r$v,sumabs[1])
    # update v
    r$v <- smr(t(x),r$u,sumabs[2])
    r$d <- t(r$u) %*% x %*% r$v
    D[c] <- r$d
  }
  return(list(r=r,D=D))
}
