#' Infer a matrix of non-negative intensities in NMF
#'
#' \code{infer_intensities} estimates a non-negative matrix \code{D} that optimizes the objective function \eqn{F = ||X - C*D||^2}
#' using per-row quadratic programming.
#'
#' @param C,X Numeric matrices.
#' @param esign A character. Keep elements of matrix \code{D} non-negative ("pos") or not ("all). By default "pos".
#' @param n.cores An integer. Number of cores to use (by default 1).
#' @return Fitted matrix \code{D}.
infer_intensities <- function(C, X, esign = "pos", n.cores = 1){
  D <- t(C) %*% C
  dmat <- X %*% C
  Amat <- diag(1, nrow(D))
  bvec <- rep(0, nrow(D))
  if (esign == "all"){
    bvec <- rep(-1e+5, nrow(D))
  }
  nr <- nrow(X)
  inten <- do.call(rbind, parallel::mclapply(1:nr, function(i){
    ft <- solve.QP(D, dmat[i, ], Amat, bvec)
    ft$solution
  },mc.cores = 20))
  colnames(inten) <- paste("comp",1:ncol(C),sep="")
  return(inten)
}
