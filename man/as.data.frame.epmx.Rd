\name{as.data.frame.epmx}
\alias{as.data.frame.epmx}
\title{Convert epmx to a data...}
\usage{\method{as.data.frame}{epmx} (x, subset, plottype="xyplot")}

\description{
  Convert epmx to a data.frame
}

\details{
  Converts an epmx to a data.frame, in the appropriate form
  for \code{\link{epmxPlotInternal}}.
}
\value{A \code{data.frame}, with a column \code{which} that denotes
the sampler, a \code{.index} column indicating the original row.}
\seealso{\code{\link{epmx}}}
\arguments{
  \item{x}{An \code{epmx} object.}
  \item{subset}{A list of length of two, the first component of which
determines the heads of the edges that are displayed, and the second
determines the tails of the edges that are displayed. The default
value \code{NULL} displays all pairs.}
  \item{plottype}{Either \code{"xyplot"} or \code{"splom"}.}
}
