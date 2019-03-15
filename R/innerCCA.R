#' Sparse Canonical Covariates Analysis (CCA)
#' 
#' This functions minimises |Cxy-udv'|2 s.t. sum(abs(u))<sumabs[1], sum(abs(v))<sumabs[2], norm(u,2)=norm(v,2) = 1 with respect to u and v. sequencially by deflation and subsequently calculates canonical variates.  
#' 
#'@param Cxy a px by py crossproduct matrix between features in X and Y
#'@param X a n by px matrix 
#'@param Y a n by py matrix 
#'@param sumabs 1 by 2 vector of L1 contraint on u and v respectively
#'@param ncomp scalar - number of components
#'
#'@example
#'
#'x <- matrix(rnorm(200),10,20)
#'results <- ssvd(x, sumabs = c(2,3))
#'
#'@export
innerCCA <- function(CxyR,X,Y,ncomp,sumabs){
  if (sqrt(ncol(X))<sumabs[1] & sqrt(ncol(Y))<sumabs[2] ){
    res <- svd(CxyR,nu = ncomp, nv = ncomp)
  } else {
    #res <- pen_nipals(CxyR,sumabs)$r
    res <- ssvd(CxyR,ncomp,sumabs)
  }
  tx <- tcrossprod(X,t(res$u))
  ty <- tcrossprod(Y,t(res$v))
  r <- diag(cor(tx,ty))
  vr <- diag(crossprod(tx,ty))
  return(list(r = r, vr = vr, tx = tx, ty = ty, Wx = res$u,Wy = res$v,d = res$d))
}
