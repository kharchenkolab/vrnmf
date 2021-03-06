% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infer_intensities.R
\name{factor_intensities}
\alias{factor_intensities}
\title{Infer a matrix of non-negative intensities in NMF with offset/nmf-offset.}
\usage{
factor_intensities(
  C,
  X,
  fit.nmf = TRUE,
  fit.factor = FALSE,
  qp.exact = FALSE,
  n.iter = 200,
  qp.iter = 10,
  rel.error.cutoff = 1e-05,
  extrapolate = TRUE,
  extrapolate.const = TRUE,
  extrapolate.convex = FALSE,
  q.factor = 1,
  verbose = TRUE,
  n.cores = 1
)
}
\arguments{
\item{C, X}{Numeric matrices.}
}
\value{
Fitted matrix \code{D}.
}
\description{
\code{factor_intensities} estimates a non-negative matrix \code{D} that optimizes the objective function \eqn{F = ||X - C*D - offset||^2},
where offset is either column-specific offset or a "1-rank nmf term": product of row vector and column vector
}
