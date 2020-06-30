# winimize W with regularized determinant
# method = det, logdet
volnmf_logdet <- function(C, X, R, R.constraint="pos", posit=FALSE, nesterov = FALSE, majorate = FALSE, beta=0.1, w.vol=1e-1, delta=1,eigen.cut = 1e-16,
                          err.cut = 1e-3, n.iter = 1e+3){

  #C <- Ctrue
  #R0 <- R <- matrix(runif(14*14,0,1),14,14)
  #R0 <- R <- Rtrue


  ## weight of volume regularization
  #w.vol <- 1e-1

  #cat('estimate volume-regularizer R.. ')

  W0 <- W <- t(R)
  FM <- solve( t(W)%*%(W) + delta*diag(1,nrow(W)) )
  H <- t(C)
  iter <- 1
  #err.cut <- 0
  obj <- vector()
  obj.violent <- 0
  err <- 1e+5
  while (err > err.cut & iter < n.iter){
    W.prev <- W
    Y <- W
    if (majorate==TRUE) {FM <- solve(t(Y)%*%Y + delta*diag(1,nrow(Y)) )}

    if (R.constraint=="pos"){
      Lip <- sqrt(sum((H%*%t(H) + w.vol*FM)^2))
      #Lip <- 1/sum((t(FM)%*%FM)^2)
      gradF <- W%*%(H%*%t(H) + w.vol*FM) - t(X)%*%t(H)
      W <- W - gradF/Lip
      W[W < 0] <- 0
    }else{
      W <- t(X)%*%t(H)%*%solve((H%*%t(H) + w.vol*FM))
    }

    if (nesterov==TRUE){
      Wfast <- W + beta*(W-W.prev)
      obj[iter] <-  sum((t(X)-W%*%H)^2) + w.vol*sum(diag(FM%*%(t(W)%*%W)))
      if (iter>=2){
        if (obj[iter] < obj[iter-1]){
          W <- Wfast
        }else{
          obj.violent <- obj.violent + 1
        }
      }
    }


    #sum((W-W.prev)^2)/sum(W^2)
    func1 <- sum((t(X) - W%*%H)^2)
    func2 <- w.vol*log(det( t(W)%*%W + delta*diag(1,nrow(Y)) ))
    func <- func1 + func2
    #print(eigen( t(Y) %*% Y + delta * diag(1, nrow(Y)) )$values)
    #print(c(func,func1,func2))
    err <- sum((W-W.prev)^2)/sum(W^2)
    #if (iter==1) { err.cut <- err }
    iter <- iter + 1
  }

  R <- t(W)
  #cat(paste('done.. ',iter,' iterations ',err,' error\n'))
  #cat(paste('violations.. ',obj.violent,'\n'))

  return(R)
}




volnmf_det <- function(C, X, R, posit=FALSE, w.vol=1e-1, eigen.cut = 1e-16,err.cut = 1e-3, n.iter = 1e+3){

  #C <- Ctrue
  #R0 <- R <- matrix(runif(14*14,0,1),14,14)
  #R0 <- R <- Rtrue


  ## weight of volume regularization
  #w.vol <- 1e-1

  #cat('estimate volume-regularizer R.. ')

  iter <- 1
  err <- 1e+5
  while (err > err.cut & iter < n.iter){
    R.prev <- R
    for (i in 1:nrow(R)){
      #for (i in 1:9){
      #print(c(iter,i))
      #i <- 1
      W <- t(R)
      WH <- W%*%t(C)

      # find quadratic approximation of determinant
      Wi <- W[,-i]
      nui <- det(t(Wi)%*%Wi)
      Xi <- t(X) - WH + as.matrix(W[,i])%*%t(as.matrix(C[,i])) #W[,i]%*%C[,i]
      #B <- diag(1,nrow(W)) - Wi%*%solve(t(Wi)%*%Wi)%*%t(Wi)
      # more efficient method to find B using null subspace analysis!
      sv <- svd( t(Wi),nu=ncol(Wi),nv=nrow(Wi) )
      nn <- sum(sv$d < eigen.cut)
      Ci <- sv$v[,(ncol(sv$v)-nn):ncol(sv$v)] # I take only last column - but how many to take??
      t(Ci)%*%Wi
      B <- Ci%*%t(Ci)
      nui*(W[,i])%*%B%*%W[,i]
      det(t(W)%*%W)

      # form quadratic form
      Qi <- sum(C[,i]^2)*diag(1,nrow(B)) + w.vol*nui*B
      #Qi <- Qi#/2
      fi <- Xi%*%C[,i]
      Amat <- diag(1,nrow(Qi));
      if (posit==TRUE ) {
        bvec <- rep(0,nrow(Qi))
      }else{
        bvec <- rep(-1e+6,nrow(Qi))
      }
      ft <- solve.QP(Dmat = Qi, dvec = fi, Amat = Amat,bvec = bvec)
      #ft$solution
      R[i,] <- ft$solution
      #print(  sum((R-Rtrue)^2)/sum(Rtrue^2))
    }
    #plot(R[i,],R0[i,])
    #sum((R-Rtrue)^2)/sum(Rtrue^2)
    err <- sum((R-R.prev)^2)/sum(R^2)
    if (iter %% 1000 == 0){print(c(iter,err))}
    iter <- iter + 1
  }

  #cat(paste('done.. ',iter,' iterations ',err,' error\n'))

  return(R)
}


