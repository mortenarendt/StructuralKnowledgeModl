#' Sparse bilinear k component model of a matrix
#' 
#' This functions minimises |x-udv'|2 s.t. sum(abs(u))<sumabs[1], sum(abs(v))<sumabs[2], norm(u,2)=norm(v,2) = 1 with respect to u and v. sequencially by deflation. 
#'@param x a n by p dataset
#'@param sumabs 1 by 2 vector ofL1 contraint on u and v respectively
#'@param niter = 100 (default) number of iterations
#'
#'@example
#'
#'x <- matrix(rnorm(200),10,20)
#'results <- ssvd(x, sumabs = c(2,3))
#'
#'@export
ssvd <- function(x,ncomp=1, sumabs = sqrt(c(nrow(x),ncol(x))), niter = 100){
  d <- v <- u <- c() 
  for (i in 1:ncomp){
    # calculate component
    r <- pen_nipals(x,sumabs,niter)
    # get 
    v <- cbind(v,r$r$v)
    u <- cbind(u,r$r$u)
    d <- c(d,r$r$d)
    # deflate
    x <- (diag(1,nrow(x)) - outer(as.vector(r$r$u),as.vector(r$r$u))) %*% x %*% (diag(1,ncol(x)) - outer(as.vector(r$r$v),as.vector(r$r$v)))
  }
  return(list(u=u,v=v,d=d))
}