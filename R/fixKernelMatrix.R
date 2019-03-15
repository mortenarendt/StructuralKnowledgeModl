#' Modify similarity matrix to a smoothing kernel
#'
#' This functions calculates scales a similarity matrix P = P' such that P 1 = 1 and P'*1 = 1
#'
#'@param P a n by n symmetric similarity matrix
#'@param thr convergency thresshold (default = 1e-5)
#'@param niter number of iterations (default = 1000)
#'
#'@example
#'
#'
#'
#'
#'@export
fixKernelMatrix <- function(P, thr = 1e-5,niter = 1000){
  conv <- 1
  n <- dim(P)[1]
  cc <- 0
  PP <- P
  while (conv>thr & cc<niter){
    cc <- cc+1
    for (i in 1:n){
      PP[i,] <- P[i,] / sum(P[i,])
      PP[,i] <- P[,i] / sum(P[,i])
    }

    conv <- norm(P - PP,'F')
    P <- PP
  }
  print(paste('Kernel Smoother updated in',cc, 'iterations at thr =',thr))
  return(PP)
}

