#' PLS calculation of weights
#'
#' Inner engine for PLS calculating weights based on a combined objective, maximising 1) covariance between X and Y, 2) description of feature to feature relationship, 3) Alternative kernel representation of sample-data and 4) with possible sparseness. 
#'@param X a (preprocessed) n by px dataset
#'@param Ys a (preprocessed) n by py dataset of responses
#'@param QQ a py by py symmetric feature kernel representing similarities between features. 
#'@param HH a px by px symmetric sample kernel representing similarities between samples. 
#'@param sumabsw a scalar contraint on the L1 norm of the weights.
#'@param niter = 20 (default) number of iterations
#'
#'@export
calc_struct_weights <- function(X,Ys,QQ,HH,aH,sumabsw = sqrt(dim(X)[2]),niter = 20){
  M <- (1-aH)*QQ %*% t(X) %*% Ys %*% t(Ys) %*% X %*% t(QQ) + aH * t(X) %*% HH %*% X
  # LS estimate
  #ws <- svd(M,nu = 1,nv = 0)$u %>% as.vector()
  # penalized estimate
  ws <- structMultMdl::ssvd(M,sumabs = c(sumabsw,sumabsw),niter = niter)$u %>% as.vector()
  # calculate eigenvalue
  d <- as.numeric(t(ws) %*% M %*% ws)
  return(list(w = ws,d = d,sumabs = sumabsw))
}
