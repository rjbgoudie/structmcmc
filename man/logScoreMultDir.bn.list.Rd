\name{logScoreMultDir.bn.list}
\alias{logScoreMultDir.bn.list}
\title{Multinomial-Dirichlet Log marginal likelihood.}
\usage{
  \method{logScoreMultDir}{bn.list} (x, data,
    hyperparameters = "bdeu", cache = new.env(hash = T),
    verbose = F, ...)
}
\arguments{
  \item{x}{An object of class "bn.list", the Bayesian
  Networks for which the marginal likelihood are computed.}

  \item{data}{A data.frame, with columns being factors
  giving the values of each random variable.}

  \item{cache}{Optionally, provide an environment with
  cached local scores for this data.}

  \item{hyperparameters}{A character vector of length one.
  Either "bdeu", "qi", "one", or "point9"}

  \item{verbose}{A logical of length 1. If true, a progress
  bar will be shown.}

  \item{...}{Further arguments (unused)}
}
\value{
  A numeric vector of length 1, giving the log marginal
  likelihood.  The environment 'cache' will also be updated
  because its scope is global.
}
\description{
  Compute the log marginal likelihood of the supplied
  Bayesian Networks.
}
\details{
  The data must be discrete. The conditional distributions
  of each random variable, conditional on its parents are
  assumed to be multinomial, with Dirichlet priors for the
  parameters.

  The notation here roughly follows Mukherjee and Speed
  (2008) Network inference using informative priors. PNAS
  105 (38) 14313-14318, doi: 10.1073/pnas.0802272105
}
\seealso{
  \code{\link{logScoreMultDir}},
  \code{\link{logScoreMultDir.bn}},
  \code{\link{logScoreMultDirOffline}},
  \code{\link{logScoreMultDirIncremental}}
}

