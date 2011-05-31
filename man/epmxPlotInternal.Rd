\name{epmxPlotInternal}
\alias{epmxPlotInternal}
\title{(Internal) Plot of cumulative edge probabilities.}
\usage{epmxPlotInternal(x, subset, plottype="xyplot")}

\description{
  (Internal) Plot of cumulative edge probabilities.
}

\details{
  Returns a xyplot/splom of the cumulative edge
  probabilities through time for bnpostmcmc.list and
  bvspostmcmc.list
}
\value{...}
\seealso{\code{\link{epmx}}}
\arguments{
  \item{x}{...}
  \item{subset}{A numeric vector. This specifies the nodes between which
the epmx should be displayed for.}
  \item{plottype}{Either \code{"xyplot"} or \code{"splom"}.}
}
