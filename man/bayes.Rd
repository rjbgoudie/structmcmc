\name{bayes}
\alias{bayes}
\title{Bayesian posterior parameter estimates.}
\usage{
  bayes(x, ...)
}
\arguments{
  \item{x}{An object}

  \item{...}{Further arguments passed to method}
}
\description{
  A generic
}
\examples{
d <- data.frame(
  a = factor(c(1, rep(3, 2), rep(1, 7))),
  b = factor(c(2, rep(1, 2), 3, 3, rep(2, 5))),
  c = factor(c(2, rep(2, 3), rep(1, 6))),
  d = factor(c(1:3, 2:3, 1, 1, 3:2, 2))
)

net <- bn(integer(0), integer(0), integer(0), c(1, 2, 3))
bayes(net, d, prior = "qi")
}
\seealso{
  \code{\link{bayes.bn}}, \code{\link{ml}}
}

