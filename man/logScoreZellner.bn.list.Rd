\name{logScoreZellner.bn.list}
\alias{logScoreZellner.bn.list}
\title{Normal Log marginal likelihood.}
\usage{\method{logScoreZellner}{bn.list} (x, data,
    cache=new.env(hash = TRUE), ...)}

\description{
  Normal Log marginal likelihood.
}

\details{
  Compute the log marginal likelihood of the supplied
  Bayesian Networks.

  The data is scored as continuous, using a form of the
  Zellner Prior.
}
\value{A numeric vector of length 1, giving the log marginal likelihood.
The environment 'cache' will also be updated because its scope is
global.}
\seealso{\code{\link{logScoreZellner.bn}}, \code{\link{logScoreZellner}}}
\arguments{
  \item{x}{An object of class "bn.list", the Bayesian Networks for which
the marginal likelihood are computed.}
  \item{data}{A matrix, with columns giving the values of each random
variable.}
  \item{cache}{Optionally, provide an environment with cached local scores
for this data.}
  \item{...}{Further arguments, currently unused.}
}
