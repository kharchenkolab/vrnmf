% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_simplex.R
\name{volnmf_simplex_col}
\alias{volnmf_simplex_col}
\title{Update of a matrix in NMF with equality contstraints on columns.}
\usage{
volnmf_simplex_col(
  X,
  R,
  C.prev = NULL,
  bound = 1,
  extrapolate = TRUE,
  err.cut = 1e-10,
  n.iter = 10000,
  qmax = 100
)
}
\arguments{
\item{X, R, C.prev}{Numeric Matrices. Matrices involved in the objective function.
Matrix \code{C.prev} serves as initialization (by default NULL).}

\item{bound}{A numeric. Equality constraint on columns of matrix \code{C} (by default 1).}

\item{extrapolate}{A boolean. Use extrapolation after local approximation (by default TRUE).}

\item{err.cut}{A numeric. Stop iterations if relative error between iterations is less than \code{err.cut} (parameter is not active now).}

\item{n.iter}{An integer. Number of iterations (by default 1,000).}
}
\value{
An updated matrix \code{C}.
}
\description{
\code{volnmf_simplex_col} finds non-negative matrix \code{C} that minimizes the objective \code{||X-C*R||^2}
under constraints that columns of C equal to 1 using local approximation with extrapolation.
}
