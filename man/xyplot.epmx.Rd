\name{xyplot.epmx}
\alias{xyplot.epmx}
\title{Plot of cumulative edge probabilities.}
\usage{
  \method{xyplot}{epmx} (x, subset = NULL)
}
\arguments{
  \item{x}{An \code{epmx} object.}

  \item{subset}{A list of length of two, the first
  component of which determines the heads of the edges that
  are displayed, and the second determines the tails of the
  edges that are displayed. The default value \code{NULL}
  displays all pairs.}
}
\value{
  An \code{xyplot}.
}
\description{
  Returns a xyplot of the cumulative edge probabilities
  through time for bnpostmcmc.list and bvspostmcmc.list
}
\seealso{
  \code{\link{epmx}}
}

