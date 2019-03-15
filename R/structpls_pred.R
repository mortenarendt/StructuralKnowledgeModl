#' Structural PLS prediction
#'
#' Predicts Yhat for a structpls model
#'@param model the output from structpls()
#'@param newdata a (preprocessed) n by px dataset of predictors
#'@param lv number of components to use from the model (default = the maximum)
#'
#'@export
structpls_pred <- function(model,newdata, lv = dim(B)[3]){
  B <- model$B
  if (lv=='all'){
    nlv <- dim(B)[3]
    Yhat <- array(NA, c(nrow(newdata),dim(model$yhat)[c(2,3)]))

    for (i in 1:nlv){
      Yhat[,,i] <- newdata %*% B[,,i]

    }
    colnames(Yhat) <- dimnames(B)[[2]]
  } else {
    # only make one-prediction from
    Yhat <- newdata %*% B[,,lv]
    colnames(Yhat) <- dimnames(B)[[2]]
  }
  return(Yhat)
}
