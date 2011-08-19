\name{mwmean.matrix}
\alias{mwmean.matrix}
\title{Moving window mean.}
\usage{
  \method{mwmean}{matrix} (x, window, ...)
}
\arguments{
  \item{x}{A matrix}

  \item{window}{...}

  \item{...}{Further arguments (unused)}
}
\value{
  A matrix of the same dimension as x, with the moving
  window values.
}
\description{
  Compute the moving window means of the columns of a
  matrix x. ie each column is treated separately
}
\seealso{
  \code{\link{mwmean}}
}

