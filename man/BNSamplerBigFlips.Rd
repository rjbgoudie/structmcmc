\name{BNSamplerBigFlips}
\alias{BNSamplerBigFlips}
\title{Big flips MCMC sampler for Bayesian networks.}
\usage{
  BNSamplerBigFlips(data, initial, prior, return =
  "network", logScoreFUN = logScoreMultDirFUN(),
  logScoreParameters = list(hyperparameters = "qi"),
  verbose = F, keepTape = F)
}
\arguments{
  \item{data}{The data.}

  \item{initial}{An object of class 'bn'. The starting
  value of the MCMC.}

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

  \item{logScoreParameters}{...}

  \item{verbose}{A logical of length 1, indicating whether
  verbose output should be printed.}

  \item{keepTape}{A logical of length 1, indicating whether
  a full log ('tape') of the MCMC sampler should be kept.
  Enabling this option can be very memory-intensive.}
}
\value{
  A function, which when called draws the next sample of
  the MCMC.
}
\description{
  Create a MCMC sampler for Bayesian Networks. The sampler
  samples Bayesian Networks (ie models).
}
\details{
  .....
}
\seealso{
  \code{\link{BNSampler}}, \code{\link{BNSamplerPT}},
  \code{\link{BNSamplerMJ}}, \code{\link{BNSamplerGrzeg}}
}

