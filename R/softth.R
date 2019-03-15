#' Soft thressholding
#' 
#' This functions truncates values in x (numerically) less than d to zero and reduces the magnitude for values larger than d by d. 
#' 
#' @param x a vector or matrix and d a scalar 
#'
#' @example 
#' softth(rnorm(100),2) 
#'
softth <- function(x, d) {return(sign(x) * pmax(0, abs(x) - d))}