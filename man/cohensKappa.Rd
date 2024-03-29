% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cohensKappa.R
\name{cohensKappa}
\alias{cohensKappa}
\title{Cohen's kappa coefficient}
\usage{
cohensKappa(x, y, adjusted = TRUE)
}
\arguments{
\item{x}{vector, coerced to factor.}

\item{y}{vector, coerced to factor.}

\item{adjusted}{logical. If \code{adjusted=FALSE}, the unadjusted accuracy
is returned instead.}
}
\value{
If adjusted=TRUE, returns the kappa coefficient [-1,1]. If adjusted=FALSE,
returns ordinary accuracy [0,1].
}
\description{
Calculates Cohen's \eqn{\kappa} (kappa) coefficient for two
equal length vectors.
}
\details{
Cohen's kappa coefficient measures the agreement between two
categorical vectors. Zero indicates that the agreement is no better than
expected by chance, while a value of one indicates perfect correspondence.
}
\note{
Cases with NA in either input vector are ignored.
}
