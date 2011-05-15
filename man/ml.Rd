\name{ml}
\alias{ml}
\title{Maximum likelihood estimates.}
\usage{ml(x, ...)}
\description{Maximum likelihood estimates.}
\details{A generic}
\seealso{\code{\link{ml.bn}}}
\arguments{\item{x}{An object}
\item{...}{Further arguments passed to method}}
\examples{d <- data.frame(
a = factor(c(1, rep(3,2), rep(1, 7))),
b = factor(c(2, rep(1, 4), rep(2, 5))),
c = factor(c(2, rep(2, 3), rep(1, 6)))
)

net <- bn(integer(0), integer(0), c(1,2))
ml(net, d)}
