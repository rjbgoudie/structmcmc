\name{exactposterior}
\alias{exactposterior}
\title{Posterior distribution on Bayesian networks.}
\usage{
  exactposterior(data,
    prior = priorUniform(empty(ncol(data), "bn")),
    logScoreFUN = logScoreMultDirFUN(),
    logScoreParameters = list(hyperparameters = "qi"),
    constraint = NULL, maxNumberParents = NULL,
    verbose = T)
}
\arguments{
  \item{data}{The data.}

  \item{prior}{A list of functions of the same length as
  \code{initial} that returns the local prior score of the
  corresponding node, given a numeric vector of parents.
  The default value \code{NULL} uses an improper uniform
  prior.}

  \item{logScoreFUN}{A list of four elements: \describe{
  \item{offline}{A function that computes the logScore of a
  Bayesian Network} \item{online}{A function that
  incrementally computes the logScore of a Bayesian
  Network} \item{local}{A function that computes the local
  logScore of a Bayesian Network} \item{prepare}{A function
  that prepares the data, and any further pre-computation
  required by the logScore functions.} }}

  \item{logScoreParameters}{A list of parameters that are
  passed to logScoreFUN.}

  \item{constraint}{A matrix of dimension ncol(data) x
  ncol(data) giving constraints to the sample space.  The
  (i, j) element is 1 if the edge i -> j is required -1 if
  the edge i -> is excluded.  0 if the edge i -> j is not
  constrained.  The diagonal of constraint must be all 0.}

  \item{maxNumberParents}{Integer of length 1. The maximum
  number of parents of any node. A \code{NULL} value uses
  the \code{ncol(data) - 1}.}

  \item{verbose}{A logical. Should a progress bar be
  displayed?}
}
\value{
  A \code{bnpost} object.
}
\description{
  Use one of a number of methods to get the posterior
  distribution
}
\seealso{
  \code{\link{posterior}}. Example priors
  \code{\link{priorGraph}}, \code{\link{priorUniform}}.
}

