# find optimum of C using ADMM
volnmf.update.C.eff <- function(X, R, C.prev = NULL, bound=1, extrapolate=TRUE, err.cut = 1e-10, n.iter = 1e+5){
  #Ctrue <- matrix( runif(192*14,0,1),nrow=192 ,ncol=14)
  #Ctrue <- apply(Ctrue,2,function(x) x/sum(x) )
  #R <- matrix( runif(ncol(C)*ncol(C),0,1),nrow=ncol(C) ,ncol=ncol(C))
  #X <- Ctrue%*%R

  #C.prev <- matrix( runif(192*14,0,1),nrow=192 ,ncol=14)
  #C.prev <- apply(C.prev,2,function(x) x/sum(x) )

  if (is.null(C.prev)){
    ft <- lm( t(X) ~ t(R) - 1) # can be estimated in closed form!
    C.prev <- t(ft$coefficients)
    C.prev[C.prev < 0] <- 0
    C.prev <- apply(C.prev,2,function(x)x/sum(x))
  }

  # precalculate matrices
  S <- R%*%t(R)
  K <- X%*%t(R)
  Lip <- sqrt(sum(S^2))

  #cat('C estimation.. ')
  err <- 1e+6
  iter <- 1
  C <- C.prev
  conv1 <- vector()
  conv2 <- vector()
  q <- c(1,(1+sqrt(5))/5)
  while( err > err.cut & iter < n.iter ){
    G <- C%*%S - K
    Chat <- C - G/Lip
    C.prev <- C
    C <- do.call(cbind,lapply(1:ncol(C),function(i){
      projection.onto.simplex(Chat[,i],bound)
    }))

    if (extrapolate==TRUE){
      extr <- (q[iter]-1)/q[iter+1]
      C <- C + extr*(C-C.prev)
    }

    #conv1[iter] <- sqrt(sum((C-Ctrue)^2))
    conv2[iter] <- sqrt(sum((C-C.prev)^2))
    #if (iter%%500==0) {print( c( iter, conv1[iter],conv2[iter] ) )}
    #if (iter%%1000==0){
    #  plot(Ctrue[,1],C[,1])
    #  legend("topleft",as.character(iter))
    #}
    iter <- iter + 1
    q[iter+1] <- (1+sqrt(1+4*q[iter]^2))/2
  }
  #cat(paste(' ',iter,' ',conv2[iter-1],'\n'))

  #plot(Ctrue[,1],C[,1])
  #apply( Ctrue-C,2,function(x) sqrt(sum(x^2)) )

  #par(mfrow=c(1,1))
  #plot(1:length(conv11),conv11,cex=0.05,pch=19,log="y")
  #points(1:length(conv1no),conv1no,cex=0.05,col="red",pch=19)
  #points(1:length(conv1),conv1,cex=0.05,col="green",pch=19)

  #par(mfrow=c(1,1))
  #plot(1:length(conv1no),conv1no,cex=0.05,pch=19,log="y")
  #points(1:length(conv1),conv1,cex=0.05,col="red",pch=19)
  #plot(1:length(conv2no),conv2no,cex=0.1,pch=19)#,log="y")
  #points(1:length(conv2),conv2,cex=0.1,pch=19,col="red")

  return(C)
}


# project unprojected vector ('unproj') onto unit probabilistic simplex
projection.onto.simplex <- function(unproj,bound){
  q <- sort(unproj,decreasing = TRUE,method="quick")
  qcum <- cumsum(q)
  mu <- (qcum-bound)/(1:length(qcum))
  cond1 <- (mu[-length(mu)] - q[-1]) > 0; cond1 <- c(cond1,TRUE)
  #cond2 <- q - mu > 0
  #ind <- which(cond1&cond2)
  ind <- which.max(cond1)
  return( pmax(0,unproj-mu[ind]) )
}



# find optimum of C using ADMM
volnmf.update.C <- function(X, R, C.prev = NULL, rho =  sum(R^2)/nrow(R), err.cut = 1e-10, n.iter = 1e+5){
  #Ctrue <- matrix( runif(192*14,0,1),nrow=192 ,ncol=14)
  #Ctrue <- apply(Ctrue,2,function(x) x/sum(x) )
  #R <- matrix( runif(ncol(C)*ncol(C),0,1),nrow=ncol(C) ,ncol=ncol(C))
  #X <- Ctrue%*%R

  if (is.null(C.prev)){
    ft <- lm( t(X) ~ t(R) - 1) # can be estimated in closed form!
    C.prev <- t(ft$coefficients)
    C.prev[C.prev < 0] <- 0
    C.prev <- apply(C.prev,2,function(x)x/sum(x))
  }

  # assign initial values
  Z <- Chat <- C <- C.prev
  V <- matrix(runif(nrow(C)*ncol(C),0,10),nrow(C),ncol(C))
  U <- runif(ncol(C),0,1)

  # precalculate inverse matrices
  Rinv <- solve(R%*%t(R) + rho*diag(1,nrow(R)))
  XTR <- X%*%t(R)
  O <- matrix(rep(1,nrow(C)^2),nrow(C),nrow(C))
  O1 <- matrix(rep(1,nrow(C)*ncol(C)),nrow(C),ncol(C))
  Oinv <- solve( O + diag(2,nrow(C)) )

  cat('C estimation.. ')
  #conv <- vector()
  #conv1 <- vector()
  err <- 1e+6
  iter <- 1
  while( err > err.cut & iter < n.iter ){
    C.prev <- C
    # update ADMM matrices
    C <- (XTR + rho*(Chat-V))%*%Rinv
    #C[C < 0] <- 0
    #Chat <- Oinv%*%( t(t(O1 - C - V) - U)  )
    Chat <- Oinv%*%( t(t(O1 + C + V + Z) - U)  )
    Z <- pmax(Chat,0)
    Chat[Chat < 0] <- 0
    V <- V + C - Chat
    U <- U + colSums(Chat)-1

    err <- sum((C-C.prev)^2)/sum(C^2)
    iter <- iter + 1
    #print(sum(Chat < 0))
    #print(c(iter, sum((C-C.prev)^2),sum((C-Ctrue)^2) ))
    #conv <- c(conv,sum((C-C.prev)^2))
    #conv1 <- c(conv1,sum((C-Ctrue)^2))
  }
  cat(paste(iter,' iterations ',err, 'error\n'))

  #plot(C[,1],Ctrue[,1],pch=19,cex=0.2);abline(a=0,b=1,col="red")
  #print(c(iter, sum((C-C.prev)^2)/sum((C)^2),sum((C-Ctrue)^2)/sum((C)^2) ))
  #iter
  return( list(C=C,Chat=Chat) )
}



# find optimum of C using ADMM
# method "qp" (default solve.QP) or "fg" (fast gradient)
volnmf.update.C1 <- function(X, R, C.prev = NULL, method="qp", meq=0, rho =  sum(R^2)/nrow(R), err.cut = 1e-10, n.iter = 1e+5){

  if (method=="qp"){
    Dmat <- (R)%*%t(R)
    Cupd <- do.call(rbind,lapply( 1:nrow(X),function(irow){
      dvec <- (R)%*%X[irow,]
      Amat <- cbind( rep(-1,nrow(R)),diag(1,nrow(R)) )
      bvec <- c(-1,rep(0,nrow(R)))
      ft <- solve.QP(Dmat,dvec,Amat,bvec,meq = meq)
      ft$solution
      #X[irow,]
      #as.numeric(t(R)%*%ft$solution  )
    }))
  }else if (method=="fg"){
    Cupd <- do.call(rbind,lapply( 1:nrow(X),function(irow){
      A <- t(R)
      Acov <- t(A)%*%A
      y <- X[irow,]
      Lip <- max(abs(eigen(A)$values^1)^2)
      Lip <- 50
      alp <- 0.5
      z <- x <- rep(1,nrow(A))/nrow(A)
      for (i in 1:1000){
        xprev <- x
        unproj <- z - as.numeric(Acov%*%z - t(A)%*%y)/Lip

        q <- sort(unproj,decreasing = TRUE)
        qcum <- cumsum(q)
        mu <- (qcum-1)/(1:length(qcum))
        cond1 <- (mu[-length(mu)] - q[-1]) > 0; cond1 <- c(cond1,TRUE)
        cond2 <- q - mu > 0
        ind <- which(cond1&cond2)
        x <- pmax(0,unproj-mu[ind])
        #alp <- (-alp^2 + sqrt(alp^4 + 4*alp^2) )/2
        beta <- alp*(1-alp)/(alp^2+alp+1)
        z <- x + beta*(x-xprev)
      }
      x
    }))
  }

  return( list(C=Cupd,Chat=Cupd) )
}


# find optimum of C
volnmf.update.Cy <- function(X, R, C.prev = NULL, n.iter = 10){
  Ctrue <- matrix( runif(192*14,0,1),nrow=192 ,ncol=14)
  Ctrue <- apply(Ctrue,2,function(x) x/sum(x) )
  R <- matrix( runif(ncol(C)*ncol(C),0,1),nrow=ncol(C) ,ncol=ncol(C))
  X0 <- Ctrue%*%R
  X <- X0 #+ rnorm(length(X0),0,0.003)

  if (is.null(C.prev)){
    ft <- lm( t(X) ~ t(R) - 1) # can be estimated in closed form!
    C.prev <- t(ft$coefficients)
    C.prev[C.prev < 0] <- 0
  }

  #C.prev <- matrix(runif(192*14,0,10),192,14)
  C <- C.prev
  C[C < 0] <- 0
  C <- apply(C,2,function(x) x/sum(x) )

  Q <- R%*%t(R)
  F <- solve(Q)
  xvec <- colSums(X)

  M <- matrix(-5,nrow=nrow(C),ncol=ncol(C))
  y <- rep(1,ncol(C))
  iter <- 1
  while(iter < 10){
    y.prev <- y
    y <- ( (R)%*%xvec - colSums(M) - Q%*%rep(1,length(y)) )/length(y)
    M <- do.call(rbind,lapply( 1:nrow(M),function(i){
      pmin(0,(R)%*%X[i,] - y)
    }))
    print(c(iter, sum((y-y.prev)^2) ))

    iter <- iter + 1
  }


  rxvec <- R%*%X[i,] <- 1
  mvec <- rep(-1,ncol(C))
  y <- rep(1,ncol(C))
  iter <- 1
  while(iter < 10){
    y.prev <- y
    y <- ( rxvec - mvec - Q%*%rep(1,length(y)) )
    M <- do.call(rbind,lapply( 1:nrow(M),function(i){
      pmin(0,(R)%*%X[i,] - y)
    }))
    print(c(iter, sum((y-y.prev)^2) ))

    iter <- iter + 1
  }



  # solution
  Csol <- do.call(rbind,lapply( 1:nrow(C),function(i){
    - as.numeric( F%*%( M[i,] + y - R%*%X[i,] )  )
  }))
  plot(Csol[,2],Ctrue[,2])

  return(C)
}



################################ old scripts and code ##########################################

# find optimum of C
volnmf.update.Cx <- function(X, R, C.prev = NULL, n.iter = 10, err.cut = 1e-3){
  #Ctrue <- matrix( runif(192*14,0,1),nrow=192 ,ncol=14)
  #Ctrue <- apply(Ctrue,2,function(x) x/sum(x) )
  #R <- matrix( runif(ncol(C)*ncol(C),0,1),nrow=ncol(C) ,ncol=ncol(C))
  #X <- Ctrue%*%R

  if (is.null(C.prev)){
    ft <- lm( t(X) ~ t(R) - 1) # can be estimated in closed form!
    C.prev <- t(ft$coefficients)
    C.prev[C.prev < 0] <- 0
  }

  #C.prev <- matrix(runif(192*14,0,10),192,14)
  C <- C.prev
  C[C < 0] <- 0
  C <- apply(C,2,function(x) x/sum(x) )

  #C <- C.prev <- Ctrue
  rm <- rowSums(R)
  rsq <- rowSums(R^2)
  iter <- 1
  err <- 1e+6
  while(iter < n.iter & err > err.cut){
    #dX <- X - C%*%R
    #dX <- matrix( runif(192*14,0,1e-10),192,14 )
    #sum(dX^2)
    #L <- ( dX%*%t(R) + t(t(C)*rsq) )
    #Q <- t(t(L)/rsq)
    #C.prev <- C


    #dX <- X - C%*%R
    #L <- ( dX%*%t(R) + t(t(C)*rsq) )
    #Q <- t(t(L)/rsq)
    #C.prev <- C

    #C <- do.call(cbind,lapply(1:ncol(C),function(k){
    for (k in 1:ncol(C)){


      dX <- X - C%*%R
      L <- ( dX%*%t(R) + t(t(C)*rsq) )
      Q <- t(t(L)/rsq)
      C.prev <- C


      Dmat <- diag(1,192)
      dvec <- Q[,k]
      Amat <- cbind( rep(1,192),diag(1,192) )
      bvec <- c(1,rep(0,192))
      ft <- solve.QP(Dmat,dvec,Amat,bvec,meq = 1)
      C[,k] <- ft$solution
    }
    #}))
    err <- sum((C-C.prev)^2)/sum(C^2)
    if(iter%%100 == 0){print(c(iter,sum((C-C.prev)^2)/sum(C^2)) )}
    #print(c(iter,sum((C-Ctrue)^2)/sum(Ctrue^2)) )
    iter <- iter+1
  }

  print(c(iter,err) )

  return(C)
}

