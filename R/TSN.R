#' Total Sum Normalization. 
#' 
#' Normalizing the rows of X to have sum == 1
#' 
#'@param X a n by p matrix 
#'
#'@export
TSN <- function(X){
  RS <- rowSums(X)
  X <- sweep(X,1,RS,'/')
  return(X)
}