#functions
# helper functions
yhat2yprob <- function(yhat){
  sumexpYhat <- apply(exp(yhat),c(1,3),sum)
  Yprob <- array(NA, dim(yhat))
  for (i in 1:dim(yhat)[2]){
    Yprob[,i,] = exp(yhat[,i,])/sumexpYhat
  }
  dimnames(Yprob) <- dimnames(yhat)
  return(Yprob)
}

simple_roc <- function(labels, scores){
  labels <- labels[order(scores, decreasing=TRUE)]
  data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
}

simple_auc <- function(TPR, FPR){
  # inputs already sorted, best scores first
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

structplsda <- function(X,y,lv,Q,H,aH,sumabsw,deflateX=FALSE){
  # lv = 3
  # Q = QQ
  # H = HH
  # aH= 0
  # sumabsw = 10
  # deflateX = T

  # make Y to dummy and scale it
  Y <- as.matrix(model.matrix(~factor(y)-1))
  Ys <- scale(Y)

  # make PLS on this system
  model <- structpls(X,Ys,lv,Q,H,aH,sumabsw,deflateX)

  model$y <- y
  # exchange yhat with probabilities
  model$yprob <- yhat2yprob(model$yhat)

  # get probability of correct class, roc, and auc for two class systems
  yprob_correct <- array(NA, dim(model$yhat)[c(1,3)])
  auc <- vector(mode = 'numeric',length = dim(model$yprob)[3])
  dfroc <- vector(mode = 'list', length = dim(model$yprob)[3])

  for (i in 1:dim(model$yprob)[3]){
    yp <- model$yprob[,,i]
    yprob_correct[,i] <- yp[Y==1]
    dfroc[[i]] <- simple_roc(y,model$yprob[,2,i])
    auc[i] <- with(dfroc[[i]], simple_auc(TPR, FPR))
  }

  model$yprob_correct <- yprob_correct
  model$auc <- auc
  model$ROC <- dfroc

  return(model)
}


mkcvindex <- function(n,k = ceiling(sqrt(n)),rep = ifelse(method =='rnd',10,1) ,method = 'rnd'){
  if (method == 'loo'){
    ID <- sample(n)
  }
  id0 <- rep(1:k,ceiling(n/k))
  id0 <- id0[1:n]
  if (method=='123') {
    ID <- id0
  }
  if (method=='111'){
    ID <- sort(id0)
  }
  if (method=='rnd'){
    ID <- matrix(NA,n,rep)
    for (i in 1:rep){
      ID[,i] <- sample(id0)
    }
  }
  return(ID)
}

check.integer <- function(x) {
  if (is.character(x)) {i <- F}
  else {
    i <- x == round(x)
  }
  return(i)
}

getsumabswgrid <- function(nsumabsw,p){
  if (nsumabsw>1 & check.integer(nsumabsw)) {nn <- nsumabsw}
  else {nn <- 0}
  sumabsw <- switch(which(c('BIC',1,nn) %in% nsumabsw),
                    -1,
                    sqrt(p),
                    seq(1,sqrt(p),length.out = nsumabsw)
  )
  return(sumabsw)
}


cv_structplsda <- function(X,y,Q,H,lv,aQ,aH,sumabsw,cvindex, deflateX=T){
  # lv <- 3
  # aQ <- 0.2
  # aH <- 0.1
  # deflateX = FALSE
  # sumabsw <- 10
  # cvindex <- CVindex[,3]
  # Q <- QQ
  # H <- HH

  # turn y into 0,1 vector
  y <- as.numeric(y)
  y <- y-min(y)
  y <- y / max(y)

  Y <- as.matrix(model.matrix(~factor(y)-1))
  ny <- length(unique(y))
  q <- aQ*Q + (1 - aQ)*diag(ncol(X))

  Yhatcv <- array(NA, c(nrow(X),ny,lv))
  mx <- max(cvindex)
  mdl_cal <- structplsda(X,y,lv,q,H,aH,sumabsw,deflateX=deflateX)
  for (i in 1:mx){
    # print(i)
    # calibration
    mdl <- structplsda(X[cvindex!=i,],y[cvindex!=i],lv,q,H[cvindex!=i,cvindex!=i],aH,sumabsw,deflateX=deflateX)
    # validatiaon
    Yhatcv[cvindex==i,,] <- structpls_pred(mdl,X[cvindex==i,],lv = 'all')
  }
  Ypropcv <- yhat2yprob(Yhatcv)

  # produce AUC metrix
  auc <- vector(mode = 'numeric',length = dim(Ypropcv)[3])
  dfroc <- vector(mode = 'list',length = dim(Ypropcv)[3])

  for (i in 1:dim(Ypropcv)[3]){
    dfroc[[i]] <- simple_roc(y,Ypropcv[,2,i])
    auc[i] <- with(dfroc[[i]], simple_auc(TPR, FPR))
  }
  results <- list()
  results$Yhatcv <- Yhatcv
  results$Ypropcv <- Ypropcv
  results$auccv <- auc
  results$dfROC <- dfroc
  results$cvindex <- cvindex
  results$calibration_model <- mdl_cal
  return(results)
}


cv_structpls <- function(X,Y,Q = diag(ncol(X)),H = diag(nrow(X)),lv = 1,aQ = 1,aH = 0,sumabsw = sqrt(ncol(X)),cvindex = mkcvindex(n = nrow(X),k = 10,rep = 1,method = 'rnd'), deflateX=FALSE){

  # X <- scale(Xmtb)
  # rownames(X) <- rownames(Xmtb)
  # Y <- scale(as.matrix(as.numeric(allfeat1$breastfeedingdays)))
  # Y <- cbind(Y,Y)
  # colnames(Y) <- paste(rep('y'),1:dim(Y)[2],sep = '')
  #
  # Q <- exp(-distQ)
  # idnan <- is.na(Q)
  # Q[idnan] <- 0
  # Q <- Q %>% doublecentering()
  # Q[idnan] <- 0
  # H <- diag(nrow(X))
  # lv <- 3
  # aQ <- 0.9
  # aH <- 0
  # sumabsw <- sqrt(ncol(X))
  # CVindex <- mkcvindex(n = nrow(X),k = 10,rep = 3,method = 'rnd')
  # cvindex <- CVindex[,2]
  # deflateX <- F



  # Q <- QQ
  # H <- HH

  ny <- dim(Y)[2]
  q <- aQ*Q + (1 - aQ)*identity(ncol(X))

  Yhatcv <- array(NA, c(nrow(X),ny,lv))
  mx <- max(cvindex)
  for (i in 1:mx){
    print(i)
    # calibration
    mdl <- structpls(X[cvindex!=i,],Y[cvindex!=i,],lv,q,H[cvindex!=i,cvindex!=i],aH,sumabsw,deflateX=deflateX)
    # validatiaon
    Yhatcv[cvindex==i,,] <- structpls_pred(mdl,X[cvindex==i,],lv = 'all')
  }

  RMSECV <- matrix(NA, ny,lv)
  for (ii in 1:lv){
    E <- Y -   Yhatcv[,,ii]
    SECV <- t(E) %*% E %>% diag
    RMSECV[,ii] <- sqrt(SECV / nrow(X))
  }

  results <- list()
  results$RMSECV <- RMSECV
  results$Yhatcv <- Yhatcv
  results$cvindex <- cvindex
  return(results)
}





## Store the svd function, but with LINPACK = T as default:
svd <- function (x, nu = min(n, p), nv = min(n, p), LINPACK = TRUE)
{
  print("LINPACK:"); print(LINPACK)  ## added so you can see it's changed
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite or missing values in 'x'")
  dx <- dim(x)
  n <- dx[1L]
  p <- dx[2L]
  if (!n || !p)
    stop("a dimension is zero")
  La.res <- La.svd(x, nu, nv)   ## your problem line
  res <- list(d = La.res$d)
  if (nu)
    res$u <- La.res$u
  if (nv) {
    if (is.complex(x))
      res$v <- Conj(t(La.res$vt))
    else res$v <- t(La.res$vt)
  }
  res
}

## Over-write current svd fn with new version:
assignInNamespace("svd", svd, "base")

