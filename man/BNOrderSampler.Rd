\name{BNOrderSampler}
\alias{BNOrderSampler}
\title{Order MCMC sampler for Bayesian Networks.}
\usage{
  BNOrderSampler(data, initial = 1:ncol(data),
    prior = priorUniform(), orderPrior = orderPriorUniform,
    return = "network", logScoreFUN = logScoreMultDirFUN(),
    logScoreParameters = list(hyperparameters = "bdeu"),
    constraint = NULL, statistics = list(nEdges = nEdges),
    maxNumberParents = NULL, moveprobs = c(0.9, 0.1, 0),
    verbose = F, keepTape = F, parentsTables = NULL,
    scoresParents = NULL)
}
\arguments{
  \item{data}{The data.}

  \item{initial}{An object of class 'bn'. The starting
  value of the MCMC.}

  \item{orderPrior}{A function that returns the prior of an
  order.}

  \item{prior}{A function that returns the prior score of
  the supplied bn.}

  \item{return}{Either "network" or "contingency".}

  \item{logScoreFUN}{A list of four elements: \describe{
  \item{offline}{A function that computes the logScore of a
  Bayesian Network} \item{online}{A function that
  incrementally computes the logScore of a Bayesian
  Network} \item{local}{A function that computes the local
  logScore of a Bayesian Network} \item{prepare}{A function
  that prepares the data, and any further pre-computation
  required by the logScore functions.} } For
  Multinomial-Dirichlet models,
  \code{\link{logScoreMultDirFUN}} returns the appropriate
  list; for Normal models with Zellner g-priors,
  \code{\link{logScoreNormalFUN}} returns the appropriate
  list.}

  \item{logScoreParameters}{A list of parameters that are
  passed to logScoreFUN.}

  \item{constraint}{A matrix of dimension ncol(data) x
  ncol(data) giving constraints to the sample space.  The
  (i, j) element is 1 if the edge i -> j is required -1 if
  the edge i -> is excluded.  0 if the edge i -> j is not
  constrained.  The diagonal of constraint must be all 0.}

  \item{statistics}{A named list of functions which should
  be applied to the current network after each step. Each
  function should accept an object of class \code{bn} and
  return a scalar output. Each item in the list must be
  named so that it can be referred to.}

  \item{maxNumberParents}{Integer of length 1. The maximum
  number of parents of any node. A \code{NULL} value gives
  the default restriction of 3.}

  \item{moveprobs}{A numeric vector of length 3. Specifies
  the probability that moves updating the parent sets of 1,
  2 and 3 nodes simultaneously.  Must sum to 1.}

  \item{verbose}{A logical of length 1, indicating whether
  verbose output should be printed.}

  \item{keepTape}{A logical of length 1, indicating whether
  a full log ('tape') of the MCMC sampler should be kept.
  Enabling this option can be very memory-intensive.}

  \item{parentsTables}{A list of tables of the form
  returned by \code{enumerateParentsTable()}}

  \item{scoresParents}{A list of the form returned by
  \code{scoreParentsTable()}}
}
\value{
  A function, which when called draws the next sample of
  the MCMC.
}
\description{
  Create a MCMC sampler for Bayesian Networks. The sampler
  samples Bayesian Networks (ie models).
}
\seealso{
  \code{\link{BNSampler}}, \code{\link{BNSamplerBigFlips}},
  \code{\link{BNSamplerPT}}, \code{\link{BNSamplerMJ}},
  \code{\link{BNSamplerGrzeg}}. Internally uses
  \code{\link{samplePair}} and \code{\link{sampleNode}}.
}

