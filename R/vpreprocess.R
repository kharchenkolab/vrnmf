#' @import Matrix
NULL

#' Preprocess the data for downstream volume analysis.
#'
#' \code{vol_preprocess} Routine normalizes the data (as requested), estimates covariance and SVD decomposition.
#'
#' @param X A numeric matrix. Covariance is estimated for column vectors of \code{X}.
#' @param col.norm A character. Specifies column normalization strategy (by default "sd"). NULL to avoid normalization.
#' @param row.norm A character. Specifies row normalization strategy (by default NULL).
#' @param pfactor A numeric A factor to normalize co-occurence matrix (by default NULL).
#' Row normalization follows column normalization. NULL to avoid normalization.
#' @return A list of objects that include normalized matrix \code{X.process}, row and column normalization factors \code{row.factors} and \code{col.factors},
#' covariance matrix \code{P0}, covariance matrix \code{P} normalized to maximum value \code{pfactor},
#' orthonormal basis \code{U} and vector of eigenvalues \code{eigens}.
#' @examples
#' small_example <- sim_factors(5, 5, 5)
#' vol <- vol_preprocess(t(small_example$X))
#'
#' @export
vol_preprocess <- function(X, col.norm = "sd", row.norm = NULL, pfactor = NULL){
  row.factors <- rep(1, nrow(X))
  col.factors <- rep(1, ncol(X))

  if (!is.null(col.norm)){
    if (col.norm == "sd"){
      col.factors <- apply(X,2,sd)
    }
  }

  if (!is.null(row.norm)) {
    if (row.norm == "sd") {
      row.factors <- apply(X,1,sd)
    }
  }

  X.process <- t(t(X) / col.factors) / row.factors
  P0 <- t(X.process) %*% X.process
  if (is.null(pfactor)) {
    pfactor <- max(P0)
  }
  P <- P0 / pfactor # how to normalize P to avoid inf problems in AchorFree?
  dimr <- svd(P)

  return( list(X.process = X.process, row.factor=row.factors, col.factors=col.factors,
               P0=P0, P = P, pfactor = pfactor, U = dimr$u, eigens = dimr$d) )
}
