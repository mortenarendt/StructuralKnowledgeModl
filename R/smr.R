#' Sparse Matrix Regression
#'
#' This functions minimises |x-uv'|2 s.t. sum(abs(u))<sumabs.
#'@param x a n by p dataset
#'@param v a 1 by p vector
#'@param a scalar value for the L1 contraint on u
#'
#'@example
#'
#'x <- matrix(rnorm(200),10,20)
#'v <- svd(x)$v[,1]
#'u <- smr(x,v,2)
#'
#'@export
smr <- function(x,v,sumabs){
  # calculate LS estimates
  u <- x %*% v
  d <- norm(as.matrix(u),'E')
  u <- u/d
  # search lambda
  if (sumabs<0){
    lambda <- searchLbic(u,x,v)
  } else {
    lambda <- searchL(u,sumabs)
  }

  u <- softth(u,lambda)
  # normalize
  u <- u/norm(u,'2')
  #print(c(norm(u,'2'),norm(u,'1'),sumabs))
  return(u)
}
