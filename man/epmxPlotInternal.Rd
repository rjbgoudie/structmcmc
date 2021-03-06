\name{epmxPlotInternal}
\alias{epmxPlotInternal}
\title{(Internal) Plot of cumulative edge probabilities.}
\usage{
  epmxPlotInternal(x, subset, plottype = "xyplot")
}
\arguments{
  \item{x}{An \code{epmx} object.}

  \item{subset}{A list of length of two, the first
  component of which determines the heads of the edges that
  are displayed, and the second determines the tails of the
  edges that are displayed. The default value \code{NULL}
  displays all pairs.}

  \item{plottype}{Either \code{"xyplot"} or
  \code{"splom"}.}
}
\value{
  An appropriate plot.
}
\description{
  Returns a xyplot/splom of the cumulative edge
  probabilities through time for bnpostmcmc.list and
  bvspostmcmc.list
}
\seealso{
  \code{\link{epmx}}
}

