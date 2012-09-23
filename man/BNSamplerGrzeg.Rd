\name{BNSamplerGrzeg}
\alias{BNSamplerGrzeg}
\title{Undocumented.}
\usage{
  BNSamplerGrzeg(data, initial, prior, return = "network",
    logScoreFUN = logScoreMultDirFUN(),
    logScoreParameters = list(hyperparameters = "qi"),
    constraint = NULL, modejumping = F, verbose = F,
    keepTape = F, grzeg = list(grzegMoveProbability = 0.5))
}
\arguments{
  \item{data}{...}

  \item{initial}{...}

  \item{prior}{...}

  \item{return}{...}

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

  \item{constraint}{...}

  \item{modejumping}{NOT IMPLEMENTED}

  \item{verbose}{...}

  \item{keepTape}{...}

  \item{grzeg}{...}
}
\description{
  method description. Not tested.
}
\seealso{
  \code{\link{BNSampler}}, \code{\link{BNSamplerPT}},
  \code{\link{BNSamplerMJ}}, \code{\link{BNGibbsSampler}}.
  Internally uses \code{\link{localPartitionFunction}}, and
  \code{\link{scoreParents}}.
}

