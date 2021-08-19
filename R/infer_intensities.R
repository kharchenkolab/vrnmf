
#' Infer a matrix of non-negative intensities in NMF
#'
#' \code{infer_intensities} estimates a non-negative matrix \code{D} that optimizes the objective function \eqn{F = ||X - C*D||^2}
#' using per-row quadratic programming.
#'
#' @param C Numeric matrices.
#' @param X Numeric matrices.
#' @param esign A character. Keep elements of matrix \code{D} non-negative ("pos") or not ("all). (default="pos")
#' @param n.cores An integer. Number of cores to use. (default=1)
#' @return Fitted matrix \code{D}.
#' 
#' @export
infer_intensities <- function(C, X, esign = "pos", n.cores = 1){
  X <- t(X)
  D <- t(C) %*% C
  dmat <- X %*% C
  Amat <- diag(1, nrow(D))
  bvec <- rep(0, nrow(D))
  if (esign == "all"){
    bvec <- rep(-1e+5, nrow(D))
  }
  nr <- nrow(X)
  inten <- do.call(cbind, parallel::mclapply(1:nr, function(i){
    ft <- quadprog::solve.QP(D, dmat[i, ], Amat, bvec)
    ft$solution
  },mc.cores = n.cores))
  rownames(inten) <- paste("comp",1:ncol(C),sep="")
  return(inten)
}



#' Infer a matrix of non-negative intensities in NMF with offset/nmf-offset.
#'
#' \code{factor_intensities} estimates a non-negative matrix \code{D} that optimizes the objective function \eqn{F = ||X - C*D - offset||^2},
#' where offset is either column-specific offset or a "1-rank nmf term": product of row vector and column vector
#' @param C Numeric matrices.
#' @param X Numeric matrices.
#' @param fit.nmf A boolean. Fit both intensities and spectrum of the offset residuals.
#' @param fit.factor A boolean. Fit only spectrum of the offset residuals (keep intensities constant across samples).
#' @param qp.exact A boolean. Estimate intensities using exact quadratic programming (qp.exact = TRUE) or inexact QP via gradient decent with extrapolation (qp.exact = FALSE).
#' @param n.iter An integer. Number of iterations.
#' @param qp.iter = 1e+1 An integer. Number of iterations of inexact QP.
#' @param rel.error.cutoff A numeric. Relative error cutoff between iterations to stop iterations.
#' @param extrapolate A boolean. Use Nesterov-like extrapolation at each iteration.
#' @param extrapolate.const A boolean. Use extrapolation scheme that adds a constant extrapolation q.factor (described below) at each iteration.
#' @param extrapolate.convex A boolean. Use Nesterov extrapolation scheme.
#' @param q.factor A numeric. Specification of a a constant extrapolation factor used in case of extrapolate.const = T.
#' @param verbose A boolean.  Print per-iteration information (by default TRUE).
#' @param n.cores An integer. Number of cores to use.
#' @return Fitted matrix \code{D}.
#' @export
factor_intensities <- function (C, X, fit.nmf = TRUE, fit.factor = FALSE, qp.exact = FALSE,
                                n.iter = 2e+2, qp.iter = 1e+1,   rel.error.cutoff = 1e-5,
                                extrapolate = TRUE, extrapolate.const = TRUE, extrapolate.convex = FALSE,
                                q.factor = 1,
                                verbose = TRUE, n.cores = 1)
{
  X <- t(X)

  ####
  #X <- as.matrix(rate.rec)
  X <- apply(X,2,function(x)x/sd(x))
  #extrapolate.const <- TRUE
  #extrapolate.convex <- FALSE
  #extrapolate.majorate <- FALSE
  #verbose = TRUE
  #qp.exact <- FALSE
  #fit.nmf <- TRUE; fit.factor <- FALSE
  #n.cores <- 1
  #n.iter <- 1e+3
  #qp.iter <- 1e+1
  #rel.error.cutoff <- 1e-5
  #q.factor <- 1
  ####

  # prepare matrices for QP
  D <- t(C) %*% C
  dmat <- X %*% C
  #if (qp.exact == TRUE){
    sv <- svd(D)
    R <- t(sv$u %*% diag(sqrt(sv$d)))
    R1 <- t(solve(R))
  #}


  iter <- 1
  inten.update <- inten <- t(matrix(1, nrow = nrow(X), ncol = ncol(C)))
  spec.offset.update <- spec.offset.old <- spec.offset <- rep(1,nrow(C))
  int.offset.update <- int.offset <- rep(1,nrow(X))
  objs <- vector(); sbio <- vector(); soff <- vector()
  ## precalculate some matrices/variables
  grad.main <- X%*%C
  X.offset.const <- grad.main%*%t(R1)
  C1 <- C%*%t(R1)
  Lip.int <- max( eigen(t(C)%*%C)$values )
  ##
  q <- c(1,(1+sqrt(5))/5)

  while (iter < n.iter){
    #print(paste("iteration.. ",iter))

    inten.old <- inten

    ### estimate intensities using exact QP
    if (qp.exact == TRUE){
      if (iter > 1) dmat <- (X - X.offset) %*% C
      #dmat1 <- dmat%*%t(R1)
      fctr <- as.numeric(spec.offset.update%*%C1)
      nr <- nrow(X)
      inten <- do.call(cbind, parallel::mclapply(1:nr, function(i) {
        ft1 <- nnls::nnls(R, X.offset.const[i,] - int.offset.update[i]*fctr)
        ft1$x
      }, mc.cores = 1, mc.preschedule = TRUE))
    }

    ### estimate intensities using inexact QP via gradient decent with extrapolation
    if (qp.exact == FALSE){
      grad.offset <- as.matrix(int.offset.update)%*%(as.matrix(spec.offset.update%*%C))
      inten.update <- inten.old
      q1 <- c(1,(1+sqrt(5))/5)
      j <- 1
      rel.error <- 1e+6
      while (j < qp.iter){
        #for (j in 1:qp.iter){
        q1[j+1] <-  (1+sqrt(1+4*q1[j]^2))/2
        inten.old1 <- inten
        grad.bio <- t(inten.update)%*%D
        grad.inten <- -grad.main + grad.offset + grad.bio
        # update inten
        inten <- inten.update - t(grad.inten) / Lip.int
        inten[inten < 0] <- 0
        # extrapolate inten
        inten.update <- inten
        extr1 <- (q1[j]-1)/q1[j+1]
        if (extrapolate == TRUE) inten.update <- inten + extr1 * (inten - inten.old1)
        rel.error <- sum( (inten - inten.old1)^2 )/sum( (inten)^2 )
        #print(paste(j,rel.error))
        j <- j + 1
      }
    }

    inten.update <- inten
    if (extrapolate == TRUE){
      extr <- (q[iter] - 1) / q[iter+1]
      inten.update <- inten + extr * (inten - inten.old)
    }


    # fit NMF or offset spectrum
    X.proc <- X - t(inten.update)%*%t(C)
    spec.offset.old <- spec.offset
    int.offset.old <- int.offset
    if (fit.nmf == TRUE){
      iter1 <- 1
      rel.error <- 1e+6
      while(iter1 < 10 & rel.error > rel.error.cutoff ){
        spec.offset.prev <- spec.offset
        spec.offset <- (t(int.offset)%*%X.proc)/sum(int.offset^2)
        spec.offset <- as.numeric(pmax(spec.offset, 0))
        int.offset <- (X.proc%*%spec.offset)/sum(spec.offset^2)
        int.offset <- as.numeric(pmax(int.offset, 0))
        rel.error <- sqrt(sum((spec.offset-spec.offset.prev)^2))/sqrt(sum((spec.offset)^2))
        iter1 <- iter1 + 1
      }
    }else if (fit.factor == TRUE){
      spec.offset <- colMeans(X.proc)
      spec.offset <- pmax(0,spec.offset)
    }

    spec.offset.update <- spec.offset
    int.offset.update <- int.offset
    if (extrapolate==TRUE){
      extr <- (q[iter]-1)/q[iter+1]
      spec.offset.update <- spec.offset + extr*(spec.offset - spec.offset.old)
      int.offset.update <- int.offset + extr*(int.offset - int.offset.old)
    }

    op <- par(no.readonly = TRUE)
    on.exit(par(op))

    if (verbose == TRUE & iter %% 5 == 0){
      X.offset <- as.matrix(int.offset)%*%t(as.matrix(spec.offset))
      X.bio <- t(inten)%*%t(C)
      X.resid <- X - X.offset - X.bio
      objs <- c(objs, sqrt(sum((X.resid)^2)))
      sbio <- c(sbio, sum(X.bio))
      soff <- c(soff, sum(X.offset))

      temppar <- par(mfrow=c(2,1),mar=c(2,2,0.5,0.5))
      on.exit(par(temppar))
      plot(1:length(objs), objs, cex = 0.5)
      plot(1:length(sbio), sbio, cex = 0.5, col = "darkred", ylim=c(min(c(sbio,soff)),max(c(sbio,soff))))
      points(1:length(soff), soff, cex = 0.5, col = "darkgreen")

      print(paste("iteration:",iter))
      print(paste("offset difference:",sqrt(sum((spec.offset-spec.offset.old)^2))/sqrt(sum((spec.offset)^2)) ))
      print( sum(X.bio) )
      print( sum(X.offset) )
      print(paste("objective:", sqrt(sum((X.resid)^2))))
    }


    iter <- iter  + 1

    if (extrapolate.const == TRUE){
      q[iter+1] <- 1 + q.factor
    }else if (extrapolate.convex == TRUE){
      q[iter+1] <- (1+sqrt(1+4*q[iter]^2))/2
    } ##else if (extrapolate.majorate == TRUE){
    ##  q[iter + 1] <- min(q.upper, q[iter + 1] * rate.q.up)
    ##}

  }

  rownames(inten) <- paste("comp", 1:ncol(C), sep = "")

  return(list(intensities = inten, spec.offset = spec.offset, int.offset = int.offset))
}

