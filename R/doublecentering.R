#' Double centring of symmetric distance/similarity matrices 
#' 
#' This functions centers a n by n symmetric matrix X resulting that coloumn and row-means are zero. #' 
#'@param M a n by n symmetric matrix 
#'@export
doublecentering <- function(M){
  R <- M*0 + rowMeans(M)
  C <- t(M*0 + colMeans(M))  # or `do.call(rbind, rep(list(colMeans(tst)),3))`
  M_double_centered = M - R - C + mean(M[])
  return(M_double_centered)  
}
