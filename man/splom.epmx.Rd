\name{splom.epmx}
\alias{splom.epmx}
\title{Plot of cumulative edge probabilities.}
\usage{\method{splom}{epmx} (x, subset)}

\description{
  Plot of cumulative edge probabilities.
}

\details{
  Returns a splom of the cumulative edge probabilities
  through time for bnpostmcmc.list and bvspostmcmc.list
}
\value{A scatter plot matrix.}
\seealso{\code{\link{epmx}}}
\arguments{
  \item{x}{An \code{epmx} object.}
  \item{subset}{A numeric vector that specifies the nodes between which
the plot will be drawn for. The default value \code{NULL} displays
all pairs.}
}
