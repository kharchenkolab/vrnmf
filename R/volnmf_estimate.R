#' Volume-regularized NMF
#'
#' \code{volnmf_main} enables volume-regularized factorization of a matrix \code{B} using the following objective function:
#' \eqn{F = ||B*Q - C*R||^2 + w.vol*volume(R)}. Matrix \code{C} is required to be non-negative and having either column or row vectors on the simplex.
#' Matrix \code{R} can optionally have non-negativity constraint. Matrix \code{Q} can optionally be identity matrix or any unitary.
#' The latter option is used to decompose co-occurence matrix \code{vol_P}.
#'
#' @param vol An output object of vol_preprocess().
#' @param B A numeric matrix. A matrix to factorize (by default NULL). If not given than matrix \code{B} is taken to be a square root decomposition of \eqn{P = B*t(B)}.
#' @param volnmf An output object of \code{volnmf.main}. An option is useful to re-estimate solution using different parameters (by default NULL).
#' @param n.comp An integer. Number of components to extract (by default 3). Defines number of columns in matrix \eqn{C}.
#' @param n.reduce An integer. Dimensional reduction of matrix B (number of columns) if taken as a square root decomposition of \code{volP} (by default equal to \code{n.comp}).
#' @param do.nmf A boolean. Estimate standard solution with \code{w.vol=0} as initialization before applying volume regularization (by default \code{TRUE}).
#' @param iter.nmf An integer. Number of iterations to get solution with \code{w.vol=0} if the former requested (by default \code{1,000}).
#' @param seed An integer. Fix seed.
#' @param domain A character. Optimize unitary rotation matrix \code{Q} ("covariance") or keep it as identity matrix (as in standard NMF). By default "covariance".
#' @param volf A character. Function that approximate volume. Can have values of "logdet" or "det" (by default "logdet").
#' @param wvol A numeric. A weight of volume-regularized term \code{volume(R)}.
#' @param delta A numeric. Logdet regularization term \code{log(det(R) + delta)} (by default 1e-8).
#' @param n.iter An integer. Number of iterations (by default \code{1,000}).
#' @param err.cut A numeric. Relative error in determinant between iterations to stop algorithm (by default \code{1e-8}).
#' @param vol.iter An integer. Number of iterations to update volume-regularized matrix \code{R} at each alternating step.
#' @param c.iter An integer. Number of iterations to update simplex matrix \code{C} at each alternating step.
#' @param extrapolate A numeric. Do Nesterov extrapolation inside blocks of R and C optimization (by default TRUE).
#' @param accelerate A numeric. Do acceleration each update after R and C blocks estimated via Nesterov-like extrapolation.
#' @param acc.C A numeric. Acceleration parameter of matrix C.
#' @param acc.R A numeric. Acceleration parameter of matrix R.
#' @param C.constraint A character. Constraint either sum of columns ("col") or sum of rows ("row) to be equal to \code{C.bound} (By default "col").
#' @param C.bound A numeric. A simplex constraint on matrix C vectors.
#' @param R.constraint A character. Set up non-negativity ("pos") constraint on elements of \code{R} (by default "pos", alternative "no").
#' @param R.majorate A boolean. Majorate logdet each iteration of \code{volnmf_logdet()} (by default FALSE).
#' @param C.init Numeric matrices. Initialization of matrices \code{C, R, Q} (by default \code{NULL}).
#' @param R.init Numeric matrices. Initialization of matrices \code{C, R, Q} (by default \code{NULL}).
#' @param Q.init Numeric matrices. Initialization of matrices \code{C, R, Q} (by default \code{NULL}).
#' @param anchor An output object of \code{AnchorFree()}. Object is used optionally to initialize matrices (by default \code{NULL}).
#' @param Ctrue A matrix. Correct matrix C if known. Useful for benchmark.
#' @param verbose A boolean. Print per-iteration information (by default FALSE).
#' @param record A numeric. Record parameters every 'record' iterations (by default \code{NULL}).
#' @param verbose.nmf A boolean. Print per-iteration information for standard NMF (by default FALSE).
#' @param record.nmf A numeric. Record parameters every 'record' iterations for standard NMF (by default \code{NULL}).
#' @param mutation.run A boolean. Assess goodness of solution using reflection test if mutation.run=TRUE (applicable only to analysis of mutation patterns).
#' @return List of objects:
#'
#' \code{C, R, Q} Factorization matrices.
#'
#' \code{C.init, R.init, Q.init} Initialization matrices for volume-regularized optimization.
#'
#' \code{C.rand, R.rand, Q.rand} Random initialization matrices for NMF optimization \code{(w.vol=0)}.
#'
#' \code{rec} a list of objects that record and store state of matrices each \code{record} iterations.
#'
#' @export
volnmf_main <- function(vol, B = NULL, volnmf = NULL, n.comp = 3, n.reduce = n.comp,
                        do.nmf=TRUE, iter.nmf = 1e+2, seed  = NULL,
                        domain = "covariance", volf = 'logdet',
                        wvol = NULL, delta = 1e-8, n.iter = 5e+2, err.cut = 1e-16,
                        vol.iter = 2e+1, c.iter = 2e+1,
                        extrapolate = TRUE, accelerate = FALSE, acc.C = 4/5, acc.R = 3/4,
                        C.constraint = "col", C.bound = 1, R.constraint = "pos", R.majorate = FALSE,
                        C.init = NULL, R.init = NULL, Q.init = NULL, anchor = NULL, Ctrue = NULL,
                        verbose = TRUE, record = 100, verbose.nmf = FALSE, record.nmf = NULL, mutation.run = FALSE){

  #B <- NULL; n.comp <- 14; n.reduce <- n.comp; volnmf <- NULL;
  #domain <- "covariance"; volf <- 'logdet';
  #wvol <- NULL; delta <- 1e-8; n.iter <- 1e+4; err.cut <- 1e-8;
  #vol.iter <- 1e+2; c.iter <- 1e+2;
  #C.constraint <- "col"; C.bound <- 1; R.constraint <- "pos";
  #C.init <- NULL; R.init <- NULL;
  #anchor = NULL;
  #frac.zero = 0.3; verbose = TRUE; record = 100
  
  if (mutation.run==FALSE) {
    rate.rec <- NULL
    xcompl <- NULL
    vol <- NULL
  }

  if (!mutation.run){
    rate.rec <- xcompl <- vol <- NULL
  }
  
  # matrix B
  if (is.null(B)){
    if (!is.null(volnmf)){
      B <- volnmf$B
    }else if(!is.null(vol)){
      B <- -vol$U[,1:n.reduce]%*%sqrt(diag(vol$eigens)[1:n.reduce,1:n.reduce])
    }
  }

  if (is.null(Q.init)) Q.init <- diag(1,nrow=n.reduce,ncol=n.comp)
  if (!is.null(seed)) set.seed(seed)

  ### initialize matrices
  if (!is.null(anchor)){   # initial matrices come from AnchorFree algorithm
    C.init <- anchor$C
    E <- anchor$E
    svE <- svd(E)
    Esq <- svE$u %*% sqrt(diag(svE$d))
    Qe <- diag(1, nrow(E))
    for (jj in 1:100){
      R <- pmax(Esq %*% Qe,0)
      Qe <- volnmf_procrustes(R, Esq)
    }
    R.init <- R
    Q.init <- Qe
    #Q.init <- diag(1, nrow = n.reduce,ncol = n.comp)
  }else{   # initial matrices are generate from uniform distribution
    if (is.null(C.init)){
      C <- matrix(runif(nrow(B)*n.comp, 0, 1), nrow = nrow(B), ncol = n.comp)

      if (C.constraint == "col"){
        C <- apply(C, 2, function(x) x / sum(x))
      }else if (C.constraint == "row"){
        C <- apply(C, 1, function(x) x / sum(x))
      }
      C.init <- C * C.bound
    }
    if (is.null(R.init)){
      R.init <- matrix(runif(n.comp*n.comp, 0, 1 / C.bound), n.comp, n.comp)

    }
  }

  C.rand <- C.init; R.rand <- R.init; Q.rand <- Q.init
  if (do.nmf == TRUE){
    message('run standard nmf.. ')
    nmf.solution <- volnmf_estimate(B, C = C.init, R = R.init, Q = Q.init,
                                    domain = domain, volf = volf,
                                    wvol = 0, delta = delta, n.iter = iter.nmf, err.cut = err.cut,
                                    #vol.iter = vol.iter / 10, c.iter = c.iter / 10,
                                    vol.iter = 5, c.iter = 5,
                                    extrapolate = extrapolate, accelerate = accelerate, acc.C = acc.C, acc.R = acc.R,
                                    C.constraint = C.constraint, C.bound = C.bound, R.constraint = R.constraint,
                                    verbose = verbose.nmf, record = record.nmf, Ctrue = Ctrue)
    message('done')
    message('\n')
    C.init <- nmf.solution$C; R.init <- nmf.solution$R; Q.init <- nmf.solution$Q
  }
  if (is.null(wvol)) wvol <- 0

  # for logdet: wvol = 0.006, for det: wvol = 5e-11 or 1e-22?
  message('run volume-regularized nmf.. ')
  vol.solution <- volnmf_estimate(B, C = C.init, R = R.init, Q = Q.init,
                                  domain = domain, volf = volf, R.majorate = R.majorate,
                                  wvol = wvol, delta = delta, n.iter = n.iter, err.cut = err.cut,
                                  vol.iter = vol.iter, c.iter = c.iter,
                                  extrapolate = extrapolate, accelerate = accelerate, acc.C = acc.C, acc.R = acc.R,
                                  C.constraint = C.constraint, C.bound = C.bound, R.constraint = R.constraint,
                                  verbose = verbose, record = record, Ctrue = Ctrue, mutation.run = mutation.run )
  message('done')
  message('\n')
  return( list( C = vol.solution$C, R = vol.solution$R, Q = vol.solution$Q,
                C.init = C.init, R.init = R.init, Q.init = Q.init,
                C.rand = C.rand, R.rand = R.rand, Q.rand = Q.rand,
                rec = vol.solution$info.record ) )
}



#' Alternating optimization of volume-regularized NMF
#'
#' \code{volnmf_estimate} provides alternating optimization of volume-regularized factorization of a matrix \code{B} using the following objective function:
#' \eqn{F = ||B*Q - C*R||^2 + w.vol*volume(R)}. Matrix \code{C} is required to be non-negative and having either column or row vectors on the simplex.
#' Matrix \code{R} can optionally have non-negativity constraint. Matrix \code{Q} can optionally be identity matrix or any unitary.
#'
#' @param B A numeric matrix. A matrix to factorize (by default NULL). If not given than matrix \code{B} is taken to be a square root decomposition of \eqn{P = B*t(B)}.
#' @param C Numeric matrices. Initial matrices for optimiztion.
#' @param R Numeric matrices. Initial matrices for optimiztion.
#' @param Q Numeric matrices. Initial matrices for optimiztion.
#' @param domain A character. Optimize unitary rotation matrix \code{Q} ("covariance") or keep it as identity matrix (as in standard NMF). By default "covariance".
#' @param volf A character. Function that approximate volume. Can have values of "logdet" or "det" (by default "logdet").
#' @param R.majorate A boolean. Majorate logdet each iteration of \code{volnmf_logdet()} (by default FALSE).
#' @param wvol A numeric. A weight of volume-regularized term \code{volume(R)}.
#' @param delta A numeric. Logdet regularization term \code{log(det(R) + delta)} (by default 1e-8).
#' @param n.iter An integer. Number of iterations (by default \code{1,000}).
#' @param err.cut A numeric. Relative error in determinant between iterations to stop algorithm (by default \code{1e-8}).
#' @param vol.iter An integer. Number of iterations to update volume-regularized matrix \code{R} at each alternating step.
#' @param c.iter An integer. Number of iterations to update simplex matrix \code{C} at each alternating step.
#' @param extrapolate A numeric. Do Nesterov extrapolation inside blocks of R and C optimization (by default TRUE).
#' @param accelerate A numeric. Do acceleration each update after R and C blocks estimated via Nesterov-like extrapolation.
#' @param acc.C A numeric. Acceleration parameter of matrix C.
#' @param acc.R A numeric. Acceleration parameter of matrix R.
#' @param C.constraint A character. Constraint either sum of columns ("col") or sum of rows ("row) to be equal to \code{C.bound} (By default "col").
#' @param C.bound A numeric. A simplex constraint on matrix C vectors.
#' @param R.constraint A character. Set up non-negativity ("pos") constraint on elements of \code{R} (by default "pos", alternative "no").
#' @param verbose A boolean. Print per-iteration information (by default FALSE)
#' @param record A numeric. Record parameters every 'record' iterations (by default \code{NULL}).
#' @param Canchor A matrix. A matrix of anchor components (unused currently). (default=NULL)
#' @param Ctrue A matrix. Correct matrix C if known. Useful for benchmark.
#' @param mutation.run A boolean. Assess goodness of solution using reflection test if mutation.run=TRUE (applicable only to analysis of mutation patterns). (default=FALSE)
#' @return List of objects:
#'
#' \code{C, R, Q}, \code{E} Factorization matrices.
#'
#' \code{iter, err} Number of iterations and relative per-iteration error \code{err} in matrix \code{C}.
#'
#' \code{info.record} a list of objects that record and store state of matrices each \code{record} iterations.
#' @export
volnmf_estimate <- function(B, C, R, Q,
                            domain = "covariance", volf = 'logdet', R.majorate = FALSE,
                            wvol = NULL, delta = 1e-8, n.iter = 1e+4, err.cut = 1e-8,
                            vol.iter = 1e+2, c.iter = 1e+2,
                            extrapolate = TRUE, accelerate = TRUE, acc.C = 4/5, acc.R = 3/4,
                            C.constraint = "col", C.bound = 1, R.constraint = "pos",
                            verbose = TRUE, record = 100, Canchor = NULL, Ctrue = NULL, mutation.run = FALSE){

  iter <- 1
  err <- 1e+5
  rvol <- vector()
  aff.mean <- vector()
  info.record <- list()
  eigens <- 1
  R.update <- R; C.update <- C
  tot.update.prev <- tot.update <- 0
  if (mutation.run==FALSE) {
    rate.rec <- NULL
    xcompl <- NULL
    vol <- NULL
  }
  
  while (iter < n.iter & err > err.cut #& (min(eigens) > 1e-16 | iter<3)
  ){
    if (domain == "covariance"){
      X <- B %*% Q
    } else{
      X <- B
    }

    ### update R
    err.prev <- sum((X - C.update %*% R)^2)
    if (volf == "logdet"){
      vol.prev <- log(det(R %*% t(R) + delta * diag(1, nrow(R))))
    } else if (volf == "det"){
      vol.prev <- det(R %*% t(R))
    }
    R.prev <- R

    if (volf == "logdet"){
      R <- volnmf_logdet(C.update, X, R.update, R.constraint = R.constraint, extrapolate = extrapolate, majorate = R.majorate,
                         w.vol = wvol, delta = delta, err.cut = 1e-100, n.iter = vol.iter)

    } else if (volf == "det"){
      R <- volnmf_det(C.update, X, R.update, posit=FALSE, w.vol=wvol, eigen.cut=1e-20, err.cut = 1e-100, n.iter = vol.iter)
    }

    err.post <- sum((X - C.update %*% R)^2)
    if (volf == "logdet"){
      vol.post <- log(det(R %*% t(R)+delta * diag(1, nrow(R))))
    } else if (volf == "det"){
      vol.post <- det(R %*% t(R))
    }
    rvol[iter] <- vol.post


    ### update C
    C.prev <- C
    if (C.constraint == "col"){
      C <- volnmf_simplex_col(X, R, C.prev = C.update, bound = C.bound, extrapolate = extrapolate,
                               err.cut = 1e-100, n.iter = c.iter)
    } else{
      C <- volnmf_simplex_row(X, R, C.prev = C.update, meq = 1)
    }
    err.post.C <- sum((X - C %*% R.update)^2)

    # accelerate C if possible
    if (accelerate == TRUE){
      C.update <- C + acc.C * (C - C.prev)
      R.update <- R + acc.R * (R - R.prev)

      C.update[C.update < 0] <- 0
      R.update[R.update < 0] <- 0

      err.update <- sum((X - C %*% R.update)^2)
      vol.update <- log(det(R.update %*% t(R.update)+delta * diag(1, nrow(R.update))))
      tot.update <- err.update + wvol*vol.update

      if (tot.update > tot.update.prev){
        C.update <- C
        R.update <- R
      }
    } else{
      C.update <- C
      R.update <- R
    }

    tot.update.prev <- tot.update
    #err.post.C.prev <- err.post.C


    ### optimize Q
    if (domain == "covariance"){
      Q <- volnmf_procrustes(C %*% R, B)
    }

    err <- sum((C - C.prev)^2) / sum(C^2)
    eigens <- eigen(R %*% t(R))$values
    aff <- 1
    if (mutation.run == TRUE){
      rownames(C) <- colnames(rate.rec)
      aff <- apply(abs(cor(C, C[xcompl, ])), 1, max)
    }else if (!is.null(Ctrue)){
      if (is.null(vol)){
        aff <- apply(cor(C, Ctrue), 1, max)
      }else{
        aff <- apply(cor(C*vol$col.factors, Ctrue), 1, max)
      }
    }
    aff.mean[iter] <- mean(aff)

    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    if (verbose == TRUE & (iter %% 100 == 0)){
      temppar <- par(mfrow=c(2,1),mar=c(4,4,1,1))
      on.exit(par(temppar))
      plot(1:iter, rvol, pch=19, cex=0.1, xlab="iteration", ylab="Vol")
      #if (!is.null(vol.ref)) {abline(h=vol.ref,col="red",lwd=1)}
      cmax <- aff.mean[length(aff.mean)]
      if (!is.null(Canchor)) {cmax <- mean(apply(abs(cor(Canchor, Canchor[xcompl, ])), 1, max))}
      plot(1:iter, aff.mean, pch=19, cex=0.1, xlab="iteration", ylab="Affinity",
           ylim=c(min(c(aff.mean, cmax)), 1))
      abline(h=cmax, col="red", lwd=1)
    }

    if (verbose==TRUE & !is.null(record)){
      if (iter %% record == 0){
        message(paste("iteration", iter, "\n"))
        message(paste("Before R update.. ","fit err:",err.prev,'vol:',wvol*vol.prev,'total:',err.prev + wvol*vol.prev,"\n" ))
        message(paste("After  R update.. ","fit err:",err.post,'vol:',wvol*vol.post,'total:',err.post + wvol*vol.post,"\n" ))
        #cat(paste("Fraction R>0: ", sum(R > -1e-10)/length(R),"\n"))
        message(paste("After  C update.. ","fit err:",err.post.C,'vol:',wvol*vol.post,'total:',err.post.C + wvol*vol.post,"\n" ))
        message(paste("Mean affinity:",mean(aff),"\n"))
        message("Affinities: ")
        message('\n')
        message(aff)
        message('\n')
        message("Eigenvalues of R%*%t(R):")
        message('\n')
        message(eigens)
        message("\n")
      }
    }

    if (!is.null(record)){
      if (iter %% record == 0){
        rec <- list(C = C, R = R, Q = Q, iter = iter)
        info.record <- c(info.record, rec)
      }
    }
    iter <- iter+1
  }

  return( list(C = C, R = R, Q = Q, iter = iter, err = err, info.record = info.record) )
}


