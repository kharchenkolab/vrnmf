
#' Simulate matrices to explores \code{vrnmf}
#'
#' \code{sim_factors} simulates non-negative factorization matrices \code{C} and \code{D}
#' under a variaty of conditions to explore factorization \eqn{X = C*D + noise}.
#'
#' @param m Integers. Size of matrices. Matrix \code{C} has a size of \code{m*r}
#' and matrix \code{D} has a size of \code{r*n}.
#' @param n Integers. Size of matrices. Matrix \code{C} has a size of \code{m*r}
#' and matrix \code{D} has a size of \code{r*n}.
#' @param r Integers. Size of matrices. Matrix \code{C} has a size of \code{m*r}
#' and matrix \code{D} has a size of \code{r*n}.
#' @param simplex A character. Either columns ("col") or rows ("row") of matrix \code{C} are projected onto unit simplex. (default="col")
#' @param distr A character. Distribution to simulate matrix entries: "unif" for uniform and "exp" for exponential distributions. (default="unif")
#' @param  frac.zeros A numeric. Fraction of zeros in matrix \code{C}. It promotes sufficient scattering of matrix column/row vectors. (default=0.4)
#' @param condition A boolean. Generate more well-conditioned matrix \code{R}. (default=FALSE)
#' @param noise A numeric. Standard deviation of gaussian noise to add. (default=0e-4)
#' @return List of simulated matrices:
#'
#' \code{X.noise}, \code{X} - noisy and original matrix \code{X} to decompose.
#'
#' \code{C}, \code{D} - factorization matrices.
#' 
#' @export
sim_factors <- function(m, n, r, simplex = "col", distr = "unif", frac.zeros = 0.4,
                      condition = FALSE, noise = 0e-4){

  # sample matrices entries from a distribution
  if (distr == "unif"){
    C <- matrix(runif(m*r,0,1), nrow = m, ncol = r)
    R <- matrix(runif(r*n,0,1), nrow = r, ncol = n)
  }else if (distr == "exp"){
    C <- matrix(rexp(m*r,1), nrow = m, ncol = r)
    R <- matrix(rexp(r*n,1), nrow = r, ncol = n)
  }

  # add zeros but avoid zero rows and columns
  if (!is.null(frac.zeros)){
    C[sample(length(C), frac.zeros * length(C))] <- 0
    C <- t(apply(C, 1, function(x){
      if (max(x) < 1e-8) x[sample(length(x), 1)] <- runif(1,0,1)
      return(x)
    }))
    C <- apply(C, 2, function(x){
      if (max(x) < 1e-8) x[sample(length(x), 1)] <- runif(1,0,1)
      return(x)
    })
  }

  # project column/row vectors onto simplex
  if (simplex=="col"){
    C <- apply(C, 2, function(x) x/sum(x))
  }else if (simplex=="row"){
    C <- t(apply(C, 1, function(x) x/sum(x)))
  }

  # condition matrix if requested
  if (condition == TRUE){
    svR <- svd(R)
    lams <- rep(1,length(svR$d))
    R1 <- svR$u %*% diag(lams) %*% t(svR$v)
    R <- R1
    R[R < 0] <- 1e-3
  }

  X <- C %*% R
  N <- matrix(rnorm(length(X), 0, noise), nrow = nrow(X), ncol = ncol(X))
  X.noise <- X + N

  return(list(X.noise = X.noise, X = X, C = C, D = R))
}



