\name{logScoreZellner}
\alias{logScoreZellner}
\title{Normal-inverse-gamma (with g-prior) Log marginal likelihood.}
\usage{logScoreZellner(x, ...)}

\description{
  Normal-inverse-gamma (with g-prior) Log marginal
  likelihood.
}

\details{
  method description
}
\seealso{\code{\link{logScoreZellner.bn}}}
\arguments{
  \item{x}{...}
  \item{...}{Further arguments passed to method}
}
\examples{data <- cbind(c(-10, -2), c(-11, -4))
net <- bn(integer(0), 1)
logScoreZellner(net, data)}
