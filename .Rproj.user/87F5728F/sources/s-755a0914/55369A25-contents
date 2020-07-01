
### input:
# vol - preprocessed object from vol.preprocess
# B - matrix to facrotize
# n.comp - number of components to extract
# volnmf - preprocessed object storing NMF results
# domain - 'covariance' or 'primary'. Formed does BQ~CR, latter does B~CR.
# volf - function that approximate volume(R) - could be 'logdet' or 'det'
# C.init, R.init - matrices of initial values
# anchor - object from AnchorFree.
# C.constarint: 'col' - sum columns to 'C.bound', 'row' - sum rows to 'C.bound'.
# R.constraint: 'pos' - non-negative values, NULL - no constraints
# verbose - print output
# record - record parameters every 'record' iterations
volnmf.main <- function(vol, B = NULL, n.comp = 3, n.reduce = n.comp, volnmf = NULL, do.nmf=TRUE, iter.nmf = 1e+3,
                        domain = "covariance", volf = 'logdet',
                        wvol = NULL, delta = 1e-8, n.iter = 1e+4, err.cut = 1e-8,
                        vol.iter = 1e+2, c.iter = 1e+2,
                        C.constraint = "col", C.bound = 1, R.constraint = "pos",R.majorate=FALSE,
                        C.init = NULL, R.init = NULL, Q.init = NULL, anchor = NULL, frac.zero = 0.3,
                        verbose = TRUE, record = NULL ){

  #B <- NULL; n.comp <- 14; n.reduce <- n.comp; volnmf <- NULL;
  #domain <- "covariance"; volf <- 'logdet';
  #wvol <- NULL; delta <- 1e-8; n.iter <- 1e+4; err.cut <- 1e-8;
  #vol.iter <- 1e+2; c.iter <- 1e+2;
  #C.constraint <- "col"; C.bound <- 1; R.constraint <- "pos";
  #C.init <- NULL; R.init <- NULL;
  #anchor = NULL;
  #frac.zero = 0.3; verbose = TRUE; record = 100

  # matrix B
  if (is.null(B)){
    if (!is.null(volnmf)){
      B <- volnmf$B
    }else if(!is.null(vol)){
      B <- -vol$U[,1:n.reduce]%*%sqrt(diag(vol$eigens)[1:n.reduce,1:n.reduce])
    }
  }

  #if (domain=="covariance"){
  if (is.null(Q.init)){
    #Q.init <- diag(1,ncol(B))
    Q.init <- diag(1,nrow=n.reduce,ncol=n.comp)
  }
  #}

  ### initialize matrices
  if (!is.null(anchor)){   # initial matrices come from AnchorFree algorithm
    C.init <- anchor$C
    E <- anchor$E
    svE <- svd(E)
    Esq <- svE$u%*%sqrt(diag(svE$d))
    Qe <- diag(1,nrow(E))
    for (jj in 1:100){
      R <- pmax( Esq%*%Qe,0 )
      Qe <- volnmf_procrustes(R,Esq)
    }
    R.init <- R
    Q.init <- Qe
    Q.init <- diag(1,nrow=n.reduce,ncol=n.comp)
  }else{   # initial matrices are generate from uniform distribution
    if (is.null(C.init)){
      rmin <- 0; cmin <- 0
      #while ( rmin < 1e-6 | cmin < 1e-6 ){
      C <- matrix(runif(nrow(B)*n.comp,0,1 ),nrow=nrow(B),ncol=n.comp)
      #C[ sample(nrow(C)*ncol(C),nrow(C)*ncol(C)*frac.zero) ] <- 0
      rmin <- min(rowSums(C)); cmin <- min(colSums(C))
      #}
      if (C.constraint=="col"){
        C <- apply(C,2,function(x) (x)/sum(x) )
      }else if (C.constraint=="row"){
        C <- apply(C,1,function(x) (x)/sum(x) )
      }
      C.init <- C*C.bound
      #range( colSums(C.init) ); range( apply(C.init,1,max) ); range( apply(C.init,2,max) )
    }
    if (is.null(R.init)){
      R.init <- matrix(runif(n.comp*n.comp,0,1/C.bound),n.comp,n.comp)
    }
  }

  C.rand <- C.init; R.rand <- R.init; Q.rand <- Q.init
  if (do.nmf==TRUE){
    cat('run standard nmf.. ')
    nmf.solution <- volnmf.estimate(B, C = C.init, R = R.init, Q = Q.init,
                                    domain = domain, volf = volf,
                                    wvol = 0, delta = delta, n.iter = iter.nmf, err.cut = err.cut,
                                    vol.iter = vol.iter/10, c.iter = c.iter/10,
                                    C.constraint = C.constraint, C.bound = C.bound, R.constraint = R.constraint,
                                    verbose = verbose, record = 20)
    cat('done');cat('\n')
    C.init <- nmf.solution$C; R.init <- nmf.solution$R; Q.init <- nmf.solution$Q
  }
  if (is.null(wvol)){ wvol <- 0 }

  # for logdet: wvol = 0.006, for det: wvol = 5e-11 or 1e-22?
  vol.solution <- volnmf.estimate(B, C = C.init, R = R.init, Q = Q.init,
                                  domain = domain, volf = volf, R.majorate=R.majorate,
                                  wvol = wvol, delta = delta, n.iter = n.iter, err.cut = err.cut,
                                  vol.iter = vol.iter, c.iter = c.iter,
                                  C.constraint = C.constraint, C.bound = C.bound, R.constraint = R.constraint,
                                  verbose = verbose, record = record )

  return( list( C = vol.solution$C, R = vol.solution$R, Q = vol.solution$Q,
                C.init = C.init, R.init = R.init, Q.init = Q.init,
                C.rand = C.rand, R.rand = R.rand, Q.rand = Q.rand,
                rec = vol.solution$info.record ) )
}




volnmf.estimate <- function(B, C = C.init, R = R.init, Q = Q.init,
                            domain = "covariance", volf = 'logdet', R.majorate=FALSE,
                            wvol = NULL, delta = 1e-8, n.iter = 1e+4, err.cut = 1e-8,
                            vol.iter = 1e+2, c.iter = 1e+2,
                            C.constraint = "col", C.bound = 1, R.constraint = "pos",
                            verbose = TRUE, record = NULL, Canchor = NULL, Ctrue = NULL ){

  iter <- 1
  err <- 1e+5
  rvol <- vector()
  aff.mean <- vector()
  info.record <- list()
  eigens <- 1

  while (iter < n.iter & err > err.cut #& (min(eigens) > 1e-16 | iter<3)
  ){
    #print(iter)
    if (domain=="covariance"){
      X <- B%*%Q
    }else{ X <- B }

    ### update R
    err.prev <- sum((X-C%*%R)^2)
    if (volf=="logdet"){
      vol.prev <- log(det(R%*%t(R)+delta*diag(1,nrow(R))))
    }else if (volf=="det"){ vol.prev <- det(R%*%t(R)) }
    R.prev <- R

    if (volf=="logdet"){
      R <- volnmf_logdet(C, X, R.prev, R.constraint=R.constraint, nesterov=FALSE, majorate = R.majorate,
                         w.vol=wvol, delta=delta,
                         err.cut = 1e-100, n.iter = vol.iter)

    }else if (volf=="det"){
      R <- volnmf_det(C, X, R.prev, posit=FALSE, w.vol=wvol, eigen.cut=1e-20,err.cut = 1e-100, n.iter = vol.iter)
    }

    err.post <- sum((X-C%*%R)^2)
    if (volf=="logdet"){
      vol.post <- log(det(R%*%t(R)+delta*diag(1,nrow(R))))
    }else if (volf=="det"){ vol.post <- det(R%*%t(R)) }
    rvol[iter] <- vol.post

    ### update C
    C.prev <- C
    if (C.constraint=="col"){
      C <- volnmf.update.C.eff(X, R, C.prev = C.prev, bound=C.bound, extrapolate = TRUE,
                               err.cut = 1e-100, n.iter = c.iter)
      #C <- volnmf.update.Cx(X, R, C.prev = C.prev, n.iter = c.iter,err.cut = 1e-100)

    }else{
      Cl <- volnmf.update.C1(X, R, C.prev = C.prev, meq=1,
                             err.cut = 1e-20, n.iter = c.iter, rho=1e+3)#,rh=0.1)
      C <- Cl$C
    }
    err.post.C <- sum((X-C%*%R)^2)

    ### optimize Q
    if (domain=="covariance"){
      Q <- volnmf_procrustes( C%*%R,B)
    }

    err <- sum((C-C.prev)^2)/sum(C^2)
    eigens <- eigen(R%*%t(R))$values
    mutation.run <- TRUE
    aff <- 1#NA
    if (mutation.run==TRUE){
      rownames(C) <- colnames(rate.rec)
      aff <- apply(abs(cor(C,C[xcompl,])),1,max)
    }else if (!is.null(Ctrue)){
      aff <- apply(cor(C,Ctrue),1,max)
    }
    aff.mean[iter] <- mean(aff)

    if (verbose==TRUE & (iter %% 100 == 0) ){
      print(str(rvol))
      print(str(aff.mean))
      par(mfrow=c(2,1),mar=c(4,4,1,1))
      plot(1:iter,rvol,pch=19,cex=0.1,xlab="iteration",ylab="Vol")
      #if (!is.null(vol.ref)) {abline(h=vol.ref,col="red",lwd=1)}
      cmax <- aff.mean[length(aff.mean)]
      if (!is.null(Canchor)) {cmax <- mean(apply(abs(cor(Canchor,Canchor[xcompl,])),1,max))}
      plot(1:iter,aff.mean,pch=19,cex=0.1,xlab="iteration",ylab="Affinity",
           ylim=c(min(c(aff.mean,cmax)),1))
      abline(h=cmax,col="red",lwd=1)
    }

    if (verbose==TRUE & (iter %% record == 0) ){
      cat(paste("iteration",iter,"\n"))
      cat(paste("Before R update.. ","fit err:",err.prev,'vol:',wvol*vol.prev,'total:',err.prev + wvol*vol.prev,"\n" ))
      cat(paste("After  R update.. ","fit err:",err.post,'vol:',wvol*vol.post,'total:',err.post + wvol*vol.post,"\n" ))
      cat(paste("Fraction R>0: ", sum(R > -1e-10)/length(R),"\n"))
      cat(paste("After  C update.. ","fit err:",err.post.C,'vol:',wvol*vol.post,'total:',err.post.C + wvol*vol.post,"\n" ))
      cat(paste("Mean affinity:",mean(aff),"\n"))
      cat("Affinities: "); cat('\n'); cat(aff); cat('\n')
      cat("Eigenvalues of R%*%t(R):"); cat('\n'); cat(eigens); cat("\n")
    }

    if (!is.null(record)){
      if (iter%%record==0){
        rec <- list(C = C, R = R, Q = Q, iter = iter)
        info.record <- c(info.record,rec)
      }
    }
    iter <- iter+1
  }

  return( list(C = C, R = R, Q = Q, iter = iter, err = err, info.record = info.record) )
}


