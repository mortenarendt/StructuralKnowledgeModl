#' Permutation Test for Regularaized Sparse Canonical Covariates Analysis
#' 
#' This functions calculates a CCA model with L1 norm sparse penalties on the canonical loadings, as well as L2 norm penalties on the structural loadings in Y. 
#' 
#'@param X a n by px matrix 
#'@param Y a n by py matrix 
#'@param R a py by py symmetric matrix indicating strucural relationship between features in Y. 
#'@param ncomp scalar - number of components
#'@param nperm scalar = 100 (default) number of permutations to be perfomed
#'@param sumabs 1 by 2 vector of L1 contraint on the canonical loadings u and v respectively
#'
#'@example
#'
#'
#'
#'
#'@export
permCCA <- function(X,Y,R,ncomp = 1,nperm = 100,sumabs){
  n <- nrow(X)
  px <- ncol(X)
  py <- ncol(Y)
  pr <- ncol(R)
  
  YR <-  tcrossprod(Y,R)
  CxyR <- crossprod(X,YR)
  # model
  r <- innerCCA(CxyR,X,Y,ncomp,sumabs)  
  # permutation loop
  vp <- rp <- matrix(NA,nrow = nperm,ncol = ncomp)
  for (i in 1:nperm){
    id <- sample(n,replace = FALSE)
    CxyRp <- crossprod(X[id,],YR)
    rrp <- innerCCA(CxyRp,X[id,],Y,ncomp,sumabs)  
    rp[i,] <- rrp$r 
    vp[i,] <- rrp$vr
  }
  p_corr <- apply(sweep(rp,2,r$r,'-')>0,2,sum)/nperm
  p_var <- apply(sweep(vp,2,r$vr,'-')>0,2,sum)/nperm
  permRes <- list(p_corr = p_corr, p_var = p_var, perm_corr = rp, perm_var = vp)
  
  return(list(model = r, permutationRes = permRes))
}
