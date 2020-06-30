
### Implementation of AnchorFree method
AnchorFree <- function(vol,n.reduce=n.comp,n.comp = 3,init=NULL,init.type="unif.both",n.iter=30,err.cut=1e-10, verbose=FALSE){
  B <- -vol$U[,1:n.reduce]%*%sqrt(diag(vol$eigens)[1:n.reduce,1:n.comp])

  M <- init
  if (is.null(M)){
    if (init.type=="diag"){
      M <- diag( runif(n.comp,-1,1)  ,n.comp)
    }else if (init.type=="similar"){
      M <- matrix(runif(n.comp*n.comp,0.9,1.1),nrow = n.comp, ncol=n.comp)
    }else if (init.type=="unif.pos"){
      M <- matrix(runif(n.comp*n.comp,0,1),nrow = n.comp, ncol=n.comp)
    }else if (init.type=="unif.both"){
      M <- matrix(runif(n.comp*n.comp,-1,1),nrow = n.comp, ncol=n.comp)
    }else if (init.type=="normal"){
      M <- matrix(rnorm(n.comp*n.comp,0,1),nrow = n.comp, ncol=n.comp)
    }else if (init.type=="ica" | init.type=="ica.pos"){
      ic <- ica::icafast((B), nc = n.comp, fun = fun)
      icM <- t(t(ic$S) + (t(solve(t(ic$M))) %*% apply(B, 2, mean))[,1])
      sgn <- (apply(icM,2,function(x) sign(sum(sign(x)*x^2)) ))
      icM <- t(t(icM)*sgn)
      if (init.type=="ica.pos"){ icM[icM < 0] <- 0 }
      ft <- lm( icM ~ B - 1 )
      M <- ft$coefficients
    }else if (init.type=="Epos"){
      preE0 <- matrix(runif(n.comp*n.comp,0,1),nrow = n.comp, ncol=n.comp)
      E0 <- preE0%*%t(preE0)
      svdE <- svd(E0)
      Mrev <- svdE$u%*%sqrt(diag(svdE$d))
      M <- solve(Mrev)
    }
  }

  M.prev <- M
  detM.prev <- det(M)
  require(lpSolveAPI)

  iter <- 1
  err <- 1e+5
  while ( iter < n.iter & abs(err) > err.cut ){
    #f <- 1; k <- 1
    for (f in 1:n.comp){
      #f <- 8
      avec <- unlist(lapply(1:n.comp,function(k){
        if (n.comp==2) {detM <- M[-k,-f]}else{detM<-det(M[-k,-f])}
        Minor <- (-1)^(k+f)*detM
      }))
      Bconstr <- rbind(B, as.numeric(rep(1,nrow(B))%*%B))

      #x <- M[,f] # + 1e-10#/183.6
      #x <- x*1.0
      #range(Bconstr[-nrow(Bconstr),]%*%x)
      #as.numeric(rep(1,nrow(B))%*%B)%*%x
      #avec%*%x
      #avec%*%get.variables(lps1)

      lps1 <- make.lp( nrow(Bconstr), n.comp)
      lp.control(lps1,sense='max')
      for(i in 1:ncol(Bconstr)){
        set.column(lps1,i,Bconstr[,i] )
      }
      set.objfn(lps1, avec)
      set.constr.type(lps1, c(rep(">=",nrow(Bconstr)-1),"=") )
      set.rhs(lps1, c(rep(0,nrow(Bconstr)-1),1) )
      set.bounds(lps1, lower = rep(-Inf,ncol(Bconstr)), columns = 1:ncol(Bconstr))
      solve(lps1)
      #get.variables(lps1)
      #get.objective(lps1)
      #avec%*%x
      #avec%*%get.variables(lps1)

      lps2 <- make.lp( nrow(Bconstr), n.comp)
      lp.control(lps2,sense='min')
      for(i in 1:ncol(Bconstr)){
        set.column(lps2,i,Bconstr[,i] )
      }
      set.objfn(lps2, avec)
      set.constr.type(lps2, c(rep(">=",nrow(Bconstr)-1),"=") )
      set.rhs(lps2, c(rep(0,nrow(Bconstr)-1),1) )
      set.bounds(lps2, lower = rep(-Inf,ncol(Bconstr)), columns = 1:ncol(Bconstr))
      solve(lps2)
      #get.variables(lps2)
      #get.objective(lps2)

      if ( abs(get.objective(lps1)) > abs(get.objective(lps2)) ){
        M[,f] <- get.variables(lps1)
      }else{  M[,f] <- get.variables(lps2) }


      #cat( c(f,det(M)) );cat('\n')
    }

    detM <- det(M)
    err <- (detM-detM.prev)/detM

    if (verbose==TRUE){
      cat( c(iter, detM, detM.prev, err) )
      cat('\n')
    }
    M.prev <- M
    detM.prev <- detM
    iter <- iter + 1
  }

  if ( abs(detM) > 1e-30 ){
    C <- B%*%M
    #rownames(C) <- colnames(rate.new)
    Ccov <- solve(t(C)%*%C)
    E <- Ccov%*%(t(C)%*%Pclean%*%C)%*%Ccov
    Ppred <- C%*%E%*%t(C)
  }else{
    C <- NA; Ccov <- NA; E <- NA; Ppred <- NA
  }


  return( list(C = C, E = E, Pest = Ppred, M = M, detM = detM, init.type = init.type) )
}
