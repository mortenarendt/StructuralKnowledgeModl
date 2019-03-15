#' Structural PLS
#'
#' Calculates a PLS model with a combined objective, maximising 1) covariance between X and Y, 2) description of feature to feature relationship, 3) Alternative kernel representation of sample-data and 4) with possible sparseness.
#'@param X a (preprocessed) n by px dataset
#'@param Y a (preprocessed) n by py dataset of responses
#'@param lv number of components
#'@param Q a py by py symmetric feature kernel representing similarities between features.
#'@param H a px by px symmetric sample kernel representing similarities between samples.
#'@param aH
#'@param sumabsw a scalar contraint on the L1 norm of the weights.
#'@param deflateX (=FALSE default) whether to deflate X in the iterative estimation of components. Y is always deflated.
#'
#'@export
structpls <- function(X,Y,lv,Q,H,aH,sumabsw,deflateX=FALSE){

  #deflateX=FALSE
  #Y <- Ys
  #lv <- 3
  #Q <- QQ
  #Q <- diag(dim(X)[2])
  #H <- HH
  #aH <- 0.3
  #aH <- 0
  #sumabsw <- sqrt(635)

  ####
  Yorg <- Y
  Xorg <- X
  ny <- dim(Y)[2]
  thr <- 1e-7
  n <- dim(X)[1]
  p <- dim(X)[2]

  NT <- c()
  basis <- xlds <- wts   <- matrix(NA,p,lv)  #  x-block weights/loadings and basis
  yscrs <- xscrs <- matrix(NA,n,lv) # y and x-block scores
  ylds  <- matrix(NA,ny,lv)           #q  y-block loadings
  for (i in 1:lv){
    # calculate weights
    rr = calc_struct_weights(X,Y,Q,H,aH,sumabsw)
    # calc loads and scores for this component
    tt <- X %*% rr$w
    normtt <- norm(tt,'2')
    NT[i] <- normtt

    vv <- pp <-  t(X) %*% tt %>% as.vector()
    vv <- vv/norm(vv,'2')
    qq <- t(Y) %*% tt %>% as.vector()
    uu <- Y %*% qq %>% as.vector()

    # deflate y
    #ttinv <- 1/(n*tt)
    ttinv <- MASS::ginv(tt)
    Dm <- diag(n) - (tt %*% ttinv)
    Y <- Dm %*% Y

    if (deflateX){
      # deflate X
      X <- Dm %*% X
    }

    # Store results
    wts[,i]   <- rr$w           #r  x-block weights
    xscrs[,i] <- tt           #t  x-block scores
    xlds[,i]  <- pp           #p  x-block loadings
    ylds[,i]  <- qq           #q  y-block loadings
    yscrs[,i] <- uu           #u  y-block scores
    basis[,i] <- vv           #v  basis of x-loadings


  }

  rownames(yscrs) <- rownames(xscrs) <- rownames(X)
  colnames(yscrs) <- colnames(xscrs) <- colnames(wts) <- colnames(xlds) <- colnames(ylds) <- paste('LV',1:lv,sep = '')
  rownames(wts) <- rownames(xlds) <- colnames(X)

  B <- array(NA, c(p,ny,lv))
  yhat <- array(NA, c(n,ny,lv))
  # calculate regressionCoefficients
  # calculate yhat for each component
  for (i in 1:lv){
    B[,,i] <- wts[,1:i] %*% MASS::ginv(xscrs[,1:i]) %*% Yorg
    yhat[,,i] <- Xorg %*% B[,,i]
  }
  dimnames(B)[[1]] <- colnames(X)
  dimnames(B)[[2]] <- dimnames(yhat)[[2]] <- colnames(Y)
  dimnames(B)[[3]] <- dimnames(yhat)[[3]] <- paste('LV',1:lv,sep = '')
  dimnames(yhat)[[1]] <- rownames(X)

  return(list(B=B,scoresX = xscrs,weigthsX = wts, loadsX = xlds, scoresY = yscrs, loadsY = ylds, yhat = yhat))
}

