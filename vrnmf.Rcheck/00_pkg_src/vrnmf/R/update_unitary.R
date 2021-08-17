#' Procrustes algorithm estimates orthonormal transformation between two matrices.
#'
#' \code{volnmf_procrustes} finds orthonormal matrix \code{Q} that minimizes objective
#' \code{||A-B*Q||^2}
#'
#' @param A Numeric Matrices. Orthonormal transformation convert matrix \code{B} in matrix \code{A}.
#' @param B Numeric Matrices. Orthonormal transformation convert matrix \code{B} in matrix \code{A}.
#' @return An optimal orthonormal tranformation matrix \code{Q}.
#' @export
volnmf_procrustes <- function(A, B){
  sv <- svd(t(A) %*% B)
  return(sv$v %*% t(sv$u))
}
