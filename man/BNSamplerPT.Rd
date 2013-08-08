\name{BNSamplerPT}
\alias{BNSamplerPT}
\title{Undocumented.}
\usage{
  BNSamplerPT(data, initial, prior, return = "network",
    logScoreFUN = logScoreMultDirFUN(),
    logScoreParameters = list(hyperparameters = "bdeu"),
    constraint = NULL, modejumping = F, verbose = F,
    keepTape = F, temp = 1,
    cache = new.env(hash = T, size = 10000L))
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

  \item{modejumping}{...}

  \item{verbose}{...}

  \item{keepTape}{...}

  \item{temp}{...}

  \item{cache}{...}
}
\description{
  method description
}
\seealso{
  \code{\link{drawPT}}. \code{\link{BNSampler}},
  \code{\link{BNSamplerMJ}}, \code{\link{BNGibbsSampler}},
  \code{\link{BNSamplerBigFlips}},
  \code{\link{BNSamplerGrzeg}}
}

