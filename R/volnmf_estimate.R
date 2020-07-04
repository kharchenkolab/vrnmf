#' Volume-regularized NMF
#'
#' \code{volnmf_main} enables volume-regularized factorization of a matrix \code{B} using the following objective function:
#' \eqn{F = ||B*Q - C*R||^2 + w.vol*volume(R)}. Matrix \code{C} is required to be non-negative and having either column or row vectors on the simplex.
#' Matrix \code{R} can optionally have non-negativity constraint. Matrix \code{Q} can optionally be identity matrix or any unitary.
#' The latter option is used to decompose co-occurence matrix \code{vol$P}.
#'
#' @param vol An output object of \code{\link{vol_preprocess()}}.
#' @param B A numeric matrix. A matrix to factorize (by default NULL). If not given than matrix \code{B} is taken to be a square root decomposition of \eqn{vol$P = B*t(B)}.
#' @param volnmf An output object of \code{volnmf.main}. An option is useful to re-estimate solution using different parameters (by default NULL).
#' @param n.comp An integer. Number of components to extract (by default 3). Defines number of columns in matrix \eqn{C}.
#' @param n.reduce An integer. Dimensional reduction of matrix B (number of columns) if taken as a square root decomposition of \code{vol$P} (by default equal to \code{n.comp}).
#' @param do.nmf A boolean. Estimate standard solution with \code{w.vol=0} as initialization before applying volume regularization (by default \code{TRUE}).
#' @param iter.nmf An integer. Number of iterations to get solution with \code{w.vol=0} if the former requested (by default \code{1,000}).
#' @param domain A character. Optimize unitary rotation matrix \code{Q} ("covariance") or keep it as identity matrix (as in standard NMF). By default "covariance".
#' @param volf A character. Function that approximate volume. Can have values of "logdet" or "det" (by default "logdet").
#' @param wvol A numeric. A weight of volume-regularized term \code{volume(R)}.
#' @param delta A numeric. Logdet regularization term \code{log(det(R) + delta)} (by default 1e-8).
#' @param n.iter An integer. Number of iterations (by default \code{1,000}).
#' @param err.cut A numeric. Relative error in determinant between iterations to stop algorithm (by default \code{1e-8}).
#' @param vol.iter An integer. Number of iterations to update volume-regularized matrix \code{R} at each alternating step.
#' @param c.iter An integer. Number of iterations to update simplex matrix \code{C} at each alternating step.
#' @param C.constraint A character. Constraint either sum of columns ("col") or sum of rows ("row) to be equal to \code{C.bound} (By default "col").
#' @param  C.bound A numeric. A simplex constraint on matrix C vectors.
#' @param R.constraint A character. Set up non-negativity ("pos") constraint on elements of \code{R} (by default "pos", alternative "no").
#' @param  R.majorate A boolean. Majorate logdet each iteration of \code{volnmf_logdet()} (by default FALSE).
#' @param  C.init,R.init,Q.init Numeric matrices. Initialization of matrices \code{C, R, Q} (by default \code{NULL}).
#' @param  anchor An output object of \code{AnchorFree()}. Object is used optionally to initialize matrices (by default \code{NULL}).
#' @param verbose A boolean. Print per-iteration information (by default FALSE)
#' @param  record A numeric. Record parameters every 'record' iterations (by default \code{NULL}).
#' @return List of objects:
#'
#' \code{C, R, Q} Factorization matrices.
#'
#' \code{C.init, R.init, Q.init} Initialization matrices for volume-regularized optimization.
#'
#' \code{C.rand, R.rand, Q.rand} Random initialization matrices for NMF optimization \code{(w.vol=0)}.
#'
#' \code{rec} a list of objects that record and store state of matrices each \code{record} iterations.
volnmf_main <- function(vol, B = NULL, volnmf = NULL, n.comp = 3, n.reduce = n.comp,
                        do.nmf=TRUE, iter.nmf = 1e+3,
                        domain = "covariance", volf = 'logdet',
                        wvol = NULL, delta = 1e-8, n.iter = 1e+4, err.cut = 1e-8,
                        vol.iter = 1e+2, c.iter = 1e+2,
                        C.constraint = "col", C.bound = 1, R.constraint = "pos", R.majorate = FALSE,
                        C.init = NULL, R.init = NULL, Q.init = NULL, anchor = NULL, Ctrue = NULL,
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

  if (is.null(Q.init)) Q.init <- diag(1,nrow=n.reduce,ncol=n.comp)

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
    cat('run standard nmf.. ')
    nmf.solution <- volnmf_estimate(B, C = C.init, R = R.init, Q = Q.init,
                                    domain = domain, volf = volf,
                                    wvol = 0, delta = delta, n.iter = iter.nmf, err.cut = err.cut,
                                    #vol.iter = vol.iter / 10, c.iter = c.iter / 10,
                                    vol.iter = 2, c.iter = 2,
                                    C.constraint = C.constraint, C.bound = C.bound, R.constraint = R.constraint,
                                    verbose = verbose, record = 20, Ctrue = Ctrue)
    cat('done'); cat('\n')
    C.init <- nmf.solution$C; R.init <- nmf.solution$R; Q.init <- nmf.solution$Q
  }
  if (is.null(wvol)) wvol <- 0

  # for logdet: wvol = 0.006, for det: wvol = 5e-11 or 1e-22?
  vol.solution <- volnmf_estimate(B, C = C.init, R = R.init, Q = Q.init,
                                  domain = domain, volf = volf, R.majorate = R.majorate,
                                  wvol = wvol, delta = delta, n.iter = n.iter, err.cut = err.cut,
                                  vol.iter = vol.iter, c.iter = c.iter,
                                  C.constraint = C.constraint, C.bound = C.bound, R.constraint = R.constraint,
                                  verbose = verbose, record = record, Ctrue = Ctrue )

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
#' @param B A numeric matrix. A matrix to factorize (by default NULL). If not given than matrix \code{B} is taken to be a square root decomposition of \eqn{vol$P = B*t(B)}.
#' @param C,R,Q Numeric matrices. Initial matrices for optimiztion.
#' @param domain A character. Optimize unitary rotation matrix \code{Q} ("covariance") or keep it as identity matrix (as in standard NMF). By default "covariance".
#' @param volf A character. Function that approximate volume. Can have values of "logdet" or "det" (by default "logdet").
#' @param  R.majorate A boolean. Majorate logdet each iteration of \code{volnmf_logdet()} (by default FALSE).
#' @param wvol A numeric. A weight of volume-regularized term \code{volume(R)}.
#' @param delta A numeric. Logdet regularization term \code{log(det(R) + delta)} (by default 1e-8).
#' @param n.iter An integer. Number of iterations (by default \code{1,000}).
#' @param err.cut A numeric. Relative error in determinant between iterations to stop algorithm (by default \code{1e-8}).
#' @param vol.iter An integer. Number of iterations to update volume-regularized matrix \code{R} at each alternating step.
#' @param c.iter An integer. Number of iterations to update simplex matrix \code{C} at each alternating step.
#' @param C.constraint A character. Constraint either sum of columns ("col") or sum of rows ("row) to be equal to \code{C.bound} (By default "col").
#' @param  C.bound A numeric. A simplex constraint on matrix C vectors.
#' @param R.constraint A character. Set up non-negativity ("pos") constraint on elements of \code{R} (by default "pos", alternative "no").
#' @param verbose A boolean. Print per-iteration information (by default FALSE)
#' @param  record A numeric. Record parameters every 'record' iterations (by default \code{NULL}).
#' @return List of objects:
#'
#' \code{C, R, Q}, \code{E} Factorization matrices.
#'
#' \code{iter, err} Number of iterations and relative per-iteration error \code{err} in matrix \code{C}.
#'
#' \code{info.record} a list of objects that record and store state of matrices each \code{record} iterations.
volnmf_estimate <- function(B, C, R, Q,
                            domain = "covariance", volf = 'logdet', R.majorate = FALSE,
                            wvol = NULL, delta = 1e-8, n.iter = 1e+4, err.cut = 1e-8,
                            vol.iter = 1e+2, c.iter = 1e+2,
                            C.constraint = "col", C.bound = 1, R.constraint = "pos",
                            verbose = TRUE, record = NULL, Canchor = NULL, Ctrue = NULL){

  iter <- 1
  err <- 1e+5
  rvol <- vector()
  aff.mean <- vector()
  info.record <- list()
  eigens <- 1

  while (iter < n.iter & err > err.cut #& (min(eigens) > 1e-16 | iter<3)
  ){
    if (domain == "covariance"){
      X <- B %*% Q
    }else{
      X <- B
    }

    ### update R
    err.prev <- sum((X - C %*% R)^2)
    if (volf == "logdet"){
      vol.prev <- log(det(R %*% t(R) + delta * diag(1, nrow(R))))
    }else if (volf == "det"){
      vol.prev <- det(R %*% t(R))
    }
    R.prev <- R

    if (volf == "logdet"){
      R <- volnmf_logdet(C, X, R.prev, R.constraint = R.constraint, nesterov = FALSE, majorate = R.majorate,
                         w.vol = wvol, delta = delta, err.cut = 1e-100, n.iter = vol.iter)

    }else if (volf == "det"){
      R <- volnmf_det(C, X, R.prev, posit=FALSE, w.vol=wvol, eigen.cut=1e-20, err.cut = 1e-100, n.iter = vol.iter)
    }

    err.post <- sum((X - C %*% R)^2)
    if (volf == "logdet"){
      vol.post <- log(det(R %*% t(R)+delta * diag(1, nrow(R))))
    }else if (volf == "det"){
      vol.post <- det(R %*% t(R))
    }
    rvol[iter] <- vol.post

    ### update C
    C.prev <- C
    if (C.constraint == "col"){
      C <- volnmf_simplex_col(X, R, C.prev = C.prev, bound = C.bound, extrapolate = TRUE,
                               err.cut = 1e-100, n.iter = c.iter)
    }else{
      C <- volnmf_simplex_row(X, R, C.prev = C.prev, meq = 1,
                             err.cut = 1e-20, n.iter = c.iter, rho=1e+3)
    }
    err.post.C <- sum((X - C %*% R)^2)

    ### optimize Q
    if (domain == "covariance"){
      Q <- volnmf_procrustes(C %*% R, B)
    }

    err <- sum((C - C.prev)^2) / sum(C^2)
    eigens <- eigen(R %*% t(R))$values
    mutation.run <- FALSE
    aff <- 1
    if (mutation.run == TRUE){
      rownames(C) <- colnames(rate.rec)
      aff <- apply(abs(cor(C, C[xcompl, ])), 1, max)
    }else if (!is.null(Ctrue)){
      aff <- apply(cor(C, Ctrue), 1, max)
    }
    aff.mean[iter] <- mean(aff)

    if (verbose == TRUE & (iter %% 100 == 0)){
      par(mfrow=c(2,1),mar=c(4,4,1,1))
      plot(1:iter, rvol, pch=19, cex=0.1, xlab="iteration", ylab="Vol")
      #if (!is.null(vol.ref)) {abline(h=vol.ref,col="red",lwd=1)}
      cmax <- aff.mean[length(aff.mean)]
      if (!is.null(Canchor)) {cmax <- mean(apply(abs(cor(Canchor, Canchor[xcompl, ])), 1, max))}
      plot(1:iter, aff.mean, pch=19, cex=0.1, xlab="iteration", ylab="Affinity",
           ylim=c(min(c(aff.mean, cmax)), 1))
      abline(h=cmax, col="red", lwd=1)
    }

    if (verbose==TRUE & (iter %% record == 0)){
      cat(paste("iteration", iter, "\n"))
      cat(paste("Before R update.. ","fit err:",err.prev,'vol:',wvol*vol.prev,'total:',err.prev + wvol*vol.prev,"\n" ))
      cat(paste("After  R update.. ","fit err:",err.post,'vol:',wvol*vol.post,'total:',err.post + wvol*vol.post,"\n" ))
      cat(paste("Fraction R>0: ", sum(R > -1e-10)/length(R),"\n"))
      cat(paste("After  C update.. ","fit err:",err.post.C,'vol:',wvol*vol.post,'total:',err.post.C + wvol*vol.post,"\n" ))
      cat(paste("Mean affinity:",mean(aff),"\n"))
      cat("Affinities: "); cat('\n'); cat(aff); cat('\n')
      cat("Eigenvalues of R%*%t(R):"); cat('\n'); cat(eigens); cat("\n")
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

