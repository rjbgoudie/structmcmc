\name{logScoreNormal.bn}
\alias{logScoreNormal.bn}
\title{Normal-inverse-gamma (with g-prior) Log marginal likelihood.}
\usage{
  \method{logScoreNormal}{bn} (x, data,
    cache = new.env(hash = T), checkInput = T, ...)
}
\arguments{
  \item{x}{An object of class "bn". The Bayesian Network by
  for which the marginal likelihood is computed.}

  \item{data}{A matrix with columns giving the value of
  each random variable.}

  \item{cache}{Optionally, provide an environment with
  cached local scores for this data.}

  \item{checkInput}{A logical of length 1, specifying
  whether to check the inputs to the function.}

  \item{...}{Further arguments, currently unused}
}
\value{
  A numeric vector of length 1, giving the log marginal
  likelihood.  The environment 'cache' will also be updated
  because its scope is global.
}
\description{
  Compute the log marginal likelihood of the supplied
  Bayesian Network.
}
\details{
  The data is scored as continuous, using a form of the
  Normal Prior.
}
\examples{
data <- cbind(c(-10, -2), c(-11, -4))
net <- bn(integer(0), 1)
logScoreNormal(net, data)
}
\seealso{
  \code{\link{logScoreNormal}},
  \code{\link{logScoreNormal.bn.list}}
}

