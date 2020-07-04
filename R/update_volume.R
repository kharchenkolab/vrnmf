#' Update volume-regularized matrix \code{R} using logdet volume approximation.
#'
#' \code{volnmf_logdet} finds matrix \code{R} that minimizes objective
#' \code{||X-C*R||^2 + w.vol*log(det(R)+delta)}.
#'
#' @param C,X,R Numeric Matrices. Matrices involved in objective function.
#' Matrix R serves as initialization.
#' @param R.constraint A character. Set up ('pos') or not ('no') non-negative constraints on matrix \code{R} (by default 'pos').
#' @param majorate A boolean. Majorate logdet each iteration (by default FALSE).
#' @param nesterov A boolean. Use Nesterov acceleration (by default FALSE, currently is not supported).
#' @param w.vol A numeric. Volume (logdet) weight in objective function.
#' @param delta A numeric. Determinant pseudocount in objective function.
#' @param err.cut A numeric. Stop algorithm if relative erro in R between iteration is less than \code{err.cut}.
#' @param n.iter An integer. Number of iterations.
#' @return An updated matrix \code{R}.
volnmf_logdet <- function(C, X, R, R.constraint = "pos",  majorate = FALSE, nesterov = FALSE,
                          w.vol = 1e-1, delta = 1, err.cut = 1e-3, n.iter = 1e+3){

  W <- t(R)
  H <- t(C)
  FM <- solve(t(W) %*% (W) + delta * diag(1,nrow(W)))

  iter <- 1
  err <- 1e+5
  #obj <- vector(); obj.violent <- 0
  while (err > err.cut & iter < n.iter){
    W.prev <- W
    Y <- W
    if (majorate == TRUE){
      FM <- solve(t(Y)%*%Y + delta*diag(1,nrow(Y)) )
    }

    if (R.constraint == "pos"){
      Lip <- sqrt(sum((H %*% t(H) + w.vol * FM)^2))
      gradF <- W %*% (H %*% t(H) + w.vol*FM) - t(X) %*% t(H)
      W <- W - gradF/Lip
      W[W < 0] <- 0
    }else{
      W <- t(X) %*% t(H) %*% solve(H %*% t(H) + w.vol * FM)
    }

    nesterov <- FALSE
    if (nesterov==TRUE){ # does not work - fix it.
      beta <- 0.1
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

    #func1 <- sum((t(X) - W%*%H)^2)
    #func2 <- w.vol*log(det( t(W)%*%W + delta*diag(1,nrow(Y)) ))
    #func <- func1 + func2
    #print(eigen( t(Y) %*% Y + delta * diag(1, nrow(Y)) )$values)
    #print(c(func,func1,func2))

    err <- sum((W-W.prev)^2)/sum(W^2)
    iter <- iter + 1
  }
  #cat(paste('done.. ',iter,' iterations ',err,' error\n'))
  #cat(paste('violations.. ',obj.violent,'\n'))

  return(t(W))
}





#' Update volume-regularized matrix \code{R} using det volume approximation
#'
#' \code{volnmf_det} finds matrix \code{R} that minimizes objective
#' \code{||X-C*R||^2 + w.vol*det(R)}
#'
#' @param C,X,R Numeric Matrices. Matrices involved in objective function.
#' Matrix R serves as initialization.
#' @param posit A boolean. Set up (TRUE) or not (FALSE) non-negative constraints on matrix \code{R} (by default FALSE).
#' @param w.vol A numeric. Volume (det) weight in objective function.
#' @param eigen.cut A numeric.
#' @param err.cut A numeric. Stop algorithm if relative erro in R between iteration is less than \code{err.cut}.
#' @param n.iter An integer. Number of iterations.
#' @return An updated matrix \code{R}.
volnmf_det <- function(C, X, R, posit=FALSE,
                       w.vol = 1e-1, eigen.cut = 1e-16, err.cut = 1e-3, n.iter = 1e+3){

  iter <- 1
  err <- 1e+5
  while (err > err.cut & iter < n.iter){
    R.prev <- R
    for (i in 1:nrow(R)){
      W <- t(R)
      WH <- W %*% t(C)

      # find quadratic approximation of determinant
      Wi <- W[, -i]
      nui <- det(t(Wi) %*% Wi)
      Xi <- t(X) - WH + as.matrix(W[, i]) %*% t(as.matrix(C[, i]))
      #B <- diag(1,nrow(W)) - Wi%*%solve(t(Wi)%*%Wi)%*%t(Wi)
      # more efficient method to find B using null subspace analysis!
      sv <- svd(t(Wi), nu = ncol(Wi), nv = nrow(Wi))
      nn <- sum(sv$d < eigen.cut)
      Ci <- sv$v[, (ncol(sv$v) - nn):ncol(sv$v)] # I take only last column - but how many to take??
      B <- Ci %*% t(Ci)
      #nui * (W[,i]) %*% B %*% W[,i]
      #det(t(W) %*% W)

      # form quadratic form
      Qi <- sum(C[, i]^2) * diag(1, nrow(B)) + w.vol * nui * B
      fi <- Xi %*% C[,i]
      Amat <- diag(1, nrow(Qi));
      if (posit == TRUE){
        bvec <- rep(0, nrow(Qi))
      }else{
        bvec <- rep(-1e+6, nrow(Qi))
      }
      ft <- solve.QP(Dmat = Qi, dvec = fi, Amat = Amat, bvec = bvec)
      R[i,] <- ft$solution
    }

    err <- sum((R-R.prev)^2)/sum(R^2)
    iter <- iter + 1
  }

  return(R)
}

