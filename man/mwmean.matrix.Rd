\name{mwmean.matrix}
\alias{mwmean.matrix}
\title{Moving window mean.}
\usage{\method{mwmean}{matrix} (x, window, ...)}

\description{
  Moving window mean.
}

\details{
  Compute the moving window means of the columns of a
  matrix x. ie each column is treated separately
}
\value{A matrix of the same dimension as x, with the moving window values.}
\seealso{\code{\link{mwmean}}}
\arguments{
  \item{x}{A matrix}
  \item{window}{...}
  \item{...}{Further arguments (unused)}
}
