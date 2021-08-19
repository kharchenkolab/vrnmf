#' @importFrom stats cor lm rexp rnorm runif sd
#' @importFrom graphics abline par points
NULL

#'
#' Non-negative tri-factorization of co-occurence matrix using minimum volume approach.
#'
#' \code{AnchorFree} method tri-factorizes (co-occurence) matrix in a product  \eqn{P ~ C*E*t(C)} of non-negative matrices \eqn{C} and \eqn{E}
#' such that matrix \eqn{E} has mininum volume and columns of matrix \eqn{C} equal to 1.
#'
#' Implementation closely follows (Fu X \emph{et al.}, IEEE Trans Pattern Anal Mach Intell., 2019).
#'
#' @param vol An output object of vol_preprocess(). The method factorizes co-occurence matrix \code{vol$P}.
#' @param n.comp An integer. Number of components to extract (by default 3). Defines number of columns in matrix \eqn{C}. (default=3)
#' @param init A numeric matrix. Initial matrix \code{M}. (default=3)
#' @param init.type A character. A strategy to randomly initialize matrix \code{M}. (default="diag") Options are to
#'
#' 1) generate diagonal unit matrix ("diag"),
#'
#' 2) use ICA solution as initialization ("ica", "ica.pos").
#'
#' or sample entries from:
#'
#' 3) uniform distribution \code{[0,1]} ("unif.pos"),
#'
#' 4) unform distribution \code{[-1,1]},
#'
#' 5) uniform distribution \code{[0.9,1.1]} ("similar"),
#'
#' 6) normal distribution \code{N(0,1)}.
#' @param n.iter An integer. Number of iterations. (default=30)
#' @param err.cut A numeric. Relative error in determinant between iterations to stop algorithm (now is not used). (default=1e-30)
#' @param verbose A boolean. Print per-iteration information (default=FALSE)
#' @return List of objects:
#'
#' \code{C}, \code{E} Factorization matrices.
#'
#' \code{Pest} Estimate of \code{vol$P} co-occurence matrix \eqn{Pest = C*E*t(C)}.
#'
#' \code{M}, \code{detM} auxiliary matrix \code{M} and its determinant.
#'
#' \code{init.type} type of initialization of matrix \code{M} that was used.
#' @examples 
#' small_example <- sim_factors(5, 5, 5)
#' vol <- vol_preprocess(t(small_example$X))
#' vol.anchor <- AnchorFree(vol)
#'
#' @export
AnchorFree <- function(vol, n.comp = 3, init = NULL, init.type = "diag",
                       n.iter = 30, err.cut = 1e-30, verbose = FALSE){

  B <- -vol$U[,1:n.comp] %*% sqrt(diag(vol$eigens)[1:n.comp,1:n.comp])
  Pclean <- B %*% t(B)

  M <- init
  if (is.null(M)){
    if (init.type == "diag"){
      #M <- diag(runif(n.comp, -1, 1), n.comp)
      M <- diag(1, n.comp)
    }else if (init.type == "similar"){
      M <- matrix(runif(n.comp*n.comp, 0.9, 1.1), nrow = n.comp, ncol = n.comp)
    }else if (init.type == "unif.pos"){
      M <- matrix(runif(n.comp*n.comp, 0, 1), nrow = n.comp, ncol = n.comp)
    }else if (init.type == "unif.both"){
      M <- matrix(runif(n.comp*n.comp, -1, 1), nrow = n.comp, ncol = n.comp)
    }else if (init.type == "normal"){
      M <- matrix(rnorm(n.comp*n.comp, 0, 1), nrow = n.comp, ncol = n.comp)
    }else if (init.type == "ica" | init.type == "ica.pos"){
      ic <- ica::icafast(B, nc = n.comp)
      icM <- t(t(ic$S) + (t(solve(t(ic$M))) %*% apply(B, 2, mean))[,1])
      sgn <- apply(icM,2,function(x) sign(sum(sign(x) * x^2)))
      icM <- t(t(icM) * sgn)
      if (init.type == "ica.pos"){ icM[icM < 0] <- 0 }
      ft <- lm( icM ~ B - 1 )
      M <- ft$coefficients
    }else if (init.type == "Epos"){
      preE0 <- matrix(runif(n.comp*n.comp, 0, 1), nrow = n.comp, ncol = n.comp)
      E0 <- preE0 %*% t(preE0)
      svdE <- svd(E0)
      Mrev <- svdE$u %*% sqrt(diag(svdE$d))
      M <- solve(Mrev)
    }
  }

  M.prev <- M
  detM.prev <- det(M)
  #require(lpSolveAPI)

  iter <- 1
  err <- 1e+5
  while ( iter < n.iter & abs(err) > err.cut ){
    for (f in 1:n.comp){
      avec <- unlist(lapply(1:n.comp, function(k){
        if (n.comp == 2) {
          detM <- M[-k, -f]
        }else{
          detM <- det(M[-k, -f])
        }
        Minor <- (-1)^(k+f) * detM
      }))
      Bconstr <- rbind(B, as.numeric(rep(1,nrow(B)) %*% B))

      #x <- M[,f] # + 1e-10#/183.6
      #x <- x*1.0
      #range(Bconstr[-nrow(Bconstr),]%*%x)
      #as.numeric(rep(1,nrow(B))%*%B)%*%x
      #avec%*%x
      #avec%*%get.variables(lps1)

      lps1 <- lpSolveAPI::make.lp( nrow(Bconstr), n.comp)
      lpSolveAPI::lp.control(lps1, sense = 'max')
      for(i in 1:ncol(Bconstr)){
        lpSolveAPI::set.column(lps1, i, Bconstr[, i])
      }
      lpSolveAPI::set.objfn(lps1, avec)
      lpSolveAPI::set.constr.type(lps1, c(rep(">=", nrow(Bconstr) - 1), "="))
      lpSolveAPI::set.rhs(lps1, c(rep(0, nrow(Bconstr) - 1), 1))
      lpSolveAPI::set.bounds(lps1, lower = rep(-Inf, ncol(Bconstr)), columns = 1:ncol(Bconstr))
      solve(lps1)
      #get.variables(lps1)
      #get.objective(lps1)
      #avec%*%x
      #avec%*%get.variables(lps1)

      lps2 <- lpSolveAPI::make.lp(nrow(Bconstr), n.comp)
      lpSolveAPI::lp.control(lps2, sense = 'min')
      for(i in 1:ncol(Bconstr)){
        lpSolveAPI::set.column(lps2, i, Bconstr[, i])
      }
      lpSolveAPI::set.objfn(lps2, avec)
      lpSolveAPI::set.constr.type(lps2, c(rep(">=", nrow(Bconstr) - 1), "="))
      lpSolveAPI::set.rhs(lps2, c(rep(0, nrow(Bconstr) - 1), 1))
      lpSolveAPI::set.bounds(lps2, lower = rep(-Inf, ncol(Bconstr)), columns = 1:ncol(Bconstr))
      solve(lps2)
      #get.variables(lps2)
      #get.objective(lps2)

      if ( abs(lpSolveAPI::get.objective(lps1)) > abs(lpSolveAPI::get.objective(lps2)) ){
        M[, f] <- lpSolveAPI::get.variables(lps1)
      }else{
        M[, f] <- lpSolveAPI::get.variables(lps2)
      }
    }

    detM <- det(M)
    err <- (detM - detM.prev) / detM

    if (verbose == TRUE){
      message( paste("iteration:", iter, "det:", detM, "prev. det:", detM.prev, "error:", err) )
      message('\n')
    }
    M.prev <- M
    detM.prev <- detM
    iter <- iter + 1
  }

  if ( abs(detM) > 1e-30 ){
    C <- B %*% M
    Ccov <- solve(t(C) %*% C)
    E <- Ccov %*% (t(C) %*% Pclean %*% C) %*% Ccov
    Ppred <- C %*% E %*% t(C)
  } else{
    C <- NA; Ccov <- NA; E <- NA; Ppred <- NA
  }


  return( list(C = C, E = E, Pest = Ppred, M = M, detM = detM, init.type = init.type) )
}
