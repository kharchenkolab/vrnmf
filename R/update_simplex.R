#' Update of a matrix in NMF with equality contstraints on columns.
#'
#' \code{volnmf_simplex_col} finds non-negative matrix \code{C} that minimizes the objective \code{||X-C*R||^2}
#' under constraints that columns of C equal to 1 using local approximation with extrapolation.
#'
#' @param X,R,C.prev Numeric Matrices. Matrices involved in the objective function.
#' Matrix \code{C.prev} serves as initialization (by default NULL).
#' @param bound A numeric. Equality constraint on columns of matrix \code{C} (by default 1).
#' @param extrapolate A boolean. Use extrapolation after local approximation (by default TRUE).
#' @param err.cut A numeric. Stop iterations if relative error between iterations is less than \code{err.cut} (parameter is not active now).
#' @param n.iter An integer. Number of iterations (by default 1,000).
#' @return An updated matrix \code{C}.
volnmf_simplex_col <- function(X, R, C.prev = NULL, bound = 1, extrapolate = TRUE,
                               err.cut = 1e-10, n.iter = 1e+4){
  if (is.null(C.prev)){
    ft <- lm( t(X) ~ t(R) - 1) # estimate in a closed form!
    C.prev <- t(ft$coefficients)
    C.prev[C.prev < 0] <- 0
    C.prev <- apply(C.prev, 2, function(x) x / sum(x) )
  }

  # precalculate matrices
  S <- R %*% t(R)
  K <- X %*% t(R)
  Lip <- sqrt(sum(S^2))

  err <- 1e+6
  iter <- 1
  C <- C.prev
  q <- c(1,(1+sqrt(5))/5)
  while(err > err.cut & iter < n.iter){
    G <- C %*% S - K
    Chat <- C - G / Lip
    C.prev <- C
    C <- do.call(cbind, lapply(1:ncol(C), function(i){
      projection_onto_simplex(Chat[, i], bound)
    }))

    if (extrapolate == TRUE){
      extr <- (q[iter] - 1) / q[iter+1]
      C <- C + extr * (C - C.prev)
    }
    iter <- iter + 1
    q[iter+1] <- (1 + sqrt(1 + 4 * q[iter]^2))/2
  }
  return(C)
}

#' Update of a matrix in NMF with equality contstraints on rows.
#'
#' \code{volnmf_simplex_row} finds non-negative matrix \code{C} that minimizes the objective \code{||X-C*R||^2}
#' under constraints that rows of C equal to 1 using per-row quadratic programming.
#'
#' @param X,R,C.prev Numeric Matrices. Matrices involved in the objective function.
#' Matrix \code{C.prev} serves as initialization (by default NULL).
#' @param meq An integer 0 or 1. Require equality (\code{meq=1}) or inequality (\code{meq=0}) constratint on rows (by default 1).
#' @return An updated matrix \code{C}.
volnmf_simplex_row <- function(X, R, C.prev = NULL, meq = 1){

  Dmat <- R %*% t(R)
  C.upd <- do.call(rbind, lapply(1:nrow(X), function(irow){
    dvec <- R %*% X[irow, ]
    Amat <- cbind(rep(-1, nrow(R)), diag(1, nrow(R)))
    bvec <- c(-1, rep(0, nrow(R)))
    ft <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
    ft$solution
  }))

  return(C.upd)
}

#' Project vector onto a probabilistic simplex.
#'
#' \code{projection_onto_simplex} projects a vector \code{unproj} onto a probabilistic simplex of sum \code{bound}.
#'
#' @param unproj A numeric vector. An unprojected vector
#' @param bound A numeric. Sum of projected vector elements.
#' @return A projected vector.
projection_onto_simplex <- function(unproj, bound){
  q <- sort(unproj, decreasing = TRUE, method = "quick")
  qcum <- cumsum(q)
  mu <- (qcum - bound) / 1:length(qcum)
  cond1 <- (mu[-length(mu)] - q[-1]) > 0
  ind <- which.max(cond1)
  return( pmax(0, unproj - mu[ind]) )
}