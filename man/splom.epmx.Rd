\name{splom.epmx}
\alias{splom.epmx}
\title{Plot of cumulative edge probabilities.}
\usage{
  \method{splom}{epmx} (x, subset = NULL)
}
\arguments{
  \item{x}{An \code{epmx} object.}

  \item{subset}{A numeric vector that specifies the nodes
  between which the plot will be drawn for. The default
  value \code{NULL} displays all pairs.}
}
\value{
  A scatter plot matrix.
}
\description{
  Returns a splom of the cumulative edge probabilities
  through time for bnpostmcmc.list and bvspostmcmc.list
}
\seealso{
  \code{\link{epmx}}
}

