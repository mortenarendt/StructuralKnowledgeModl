#' Find parameter for Soft thressholding using Bayes Information Criterion
#'
#' This functions searches the optimal penalty based on BIC
#'@param w an initial guess
#'@param M a cross product matrix
#'@param n number of observations
#'@param m=30 (default) number of splits to tests
#'@export
searchLbic <- function(u,x,v,m = 30, n = dim(x)[1]){
  L = seq(0,max(abs(u))*0.99,length.out = m)
  BIC2 <- matrix(NA, nrow = m,ncol = 2)
  for (i in 1:m){
    #shrink by L[i]
    us <- structMultMdl:::softth(u,L[i])
    us <- us / norm(us,'2')
    # evaluate fit
    BIC2[i,1] <- -log(t(us) %*% x  %*% v / n)
    BIC2[i,2] = log(n)/n*sum(us!=0)
  }
  BIC = rowSums(BIC2)
  l <- L[which.min(BIC)]
  return(l)
}

# searchLbic <- function(w,M,n, m = 30){
#   L = seq(0,max(abs(w))*0.99,length.out = m)
#   BIC2 <- matrix(NA, nrow = m,ncol = 2)
#   for (i in 1:m){
#     # shrink the vector by L(i)
#     ws <- structMultMdl:::softth(w,L[i])
#     ws <- ws / norm(ws,'2')
#     BIC2[i,1] <- -log(t(ws) %*% M  %*% ws / n)
#     BIC2[i,2] = log(n)/n*sum(ws!=0)
#   }
#   BIC = rowSums(BIC2)
#   l <- L[which.min(BIC)]
#   return(l)
# }
