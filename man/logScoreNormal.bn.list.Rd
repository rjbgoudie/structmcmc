\name{logScoreNormal.bn.list}
\alias{logScoreNormal.bn.list}
\title{Normal-inverse-gamma (with g-prior) Log marginal likelihood.}
\usage{
  \method{logScoreNormal}{bn.list} (x, data,
    cache = new.env(hash = TRUE), ...)
}
\arguments{
  \item{x}{An object of class "bn.list", the Bayesian
  Networks for which the marginal likelihood are computed.}

  \item{data}{A matrix, with columns giving the values of
  each random variable.}

  \item{cache}{Optionally, provide an environment with
  cached local scores for this data.}

  \item{...}{Further arguments, currently unused.}
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
  The data is scored as continuous, using a form of the
  Normal Prior.
}
\seealso{
  \code{\link{logScoreNormal.bn}},
  \code{\link{logScoreNormal}}
}

