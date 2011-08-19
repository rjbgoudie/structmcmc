\name{ml.bn}
\alias{ml.bn}
\title{Maximum likelihood estimates for parameters of a BN.}
\usage{
  \method{ml}{bn} (x, data, nodes = seq_along(x),
  regularisation = NaN, cache = new.env(hash = T), ...)
}
\arguments{
  \item{x}{The Bayesian Network. An object of class 'bn'}

  \item{data}{A data frame}

  \item{nodes}{A subset of 1, ..., \code{nNodes(x)}. A
  numeric vector.}

  \item{regularisation}{One of \code{NaN}, \code{"qi"}, or
  a numeric vector of length 1. Supplying \code{NaN} will
  place a \code{NaN} in any parts of the table where there
  is no information. Supplying \code{"qi"} will add a
  factor of \code{1/(number of parents of node)} to each
  value (as in Bayesian inference). Supplying a number will
  add this to each value.}

  \item{cache}{A cache}

  \item{...}{Further arguments (unused)}
}
\value{
  A list of length \code{nNodes(x)}. Each component is a
  list containing components for each configuration of that
  node's parents in the Bayesian Network \code{x}. Each of
  these components is a numeric vector of probabilities
  that sum to 1, labelled with the levels of the relevant
  node.
}
\description{
  Computes the maximum likelihood estimates for the
  parameters of a Bayesian Network. These are just the
  proportions of each category for each configuration of
  the parents of a node in the Bayesian Network.
}
\details{
  If, for a particular configuration of the parents, we
  have no observations, then by default, \code{NaN} is
  returned.
}
\seealso{
  \code{\link{ml}}
}

