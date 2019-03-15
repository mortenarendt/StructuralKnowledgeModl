#' Regularaized Sparse Canonical Covariates Analysis
#'
#' This functions calculates a CCA model with L1 norm sparse penalties on the canonical loadings, as well as L2 norm penalties on the structural loadings in Y.
#'
#'@param X a n by px matrix
#'@param Y a n by py matrix
#'@param R a py by py symmetric matrix indicating strucural relationship between features in Y.
#'@param ncomp = 1 (default) - number of components
#'@param nperm = 100 (default) - number of permutations to be perfomed
#'@param sumabs  = sqrt(c(px,py)) - vector of L1 contraint on the canonical loadings u and v respectively
#'@param na = 30 (default) - indicating the grid of structural penalties based on R to tests or a vector of specific penalties.
#'
#'@example
#'
#'
#'
#'
#'@export
pCCA <- function(X,Y,R = diag(1,ncol(Y)),ncomp=1,sumabs = c(sqrt(ncol(X)),sqrt(ncol(Y))),na=30,nperm=100){
  pr <- ncol(Y)
  if (na>1){A <- seq(0,1,length.out = na)}
  else {A <- na; na <- 1}

  p_corr <- p_var <- matrix(NA,nrow = na,ncol = ncomp)
  results <- list()
  c <- 0
  for (a in A){
    c <- c+1
    Ra <- R*a + diag(1-a,pr)
    #YR <-  tcrossprod(Y,Ra)
    #Cxy <- crossprod(X,Y)
    #CxyR <- tcrossprod(Cxy,Ra)
    #CxyR <- crossprod(X,YR)
    results[[c]] <- permCCA(X,Y,Ra,ncomp,nperm,sumabs)
    p_corr[c,] <- results[[c]]$permutationRes$p_corr
    p_var[c,] <- results[[c]]$permutationRes$p_var
  }

  colnames(p_corr) <- paste('p_corr_cmp',1:ncomp,sep = '')
  colnames(p_var) <- paste('p_var_cmp',1:ncomp,sep = '')
  res_df <- data.frame(A = A,p_corr,p_var)
  return(list(res_df = res_df, a_results = results))
}


