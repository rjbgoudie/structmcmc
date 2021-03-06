\name{levelplot.ep.list}
\alias{levelplot.ep.list}
\title{Levelplot of posterior edge probabilities.}
\usage{
  \method{levelplot}{ep.list} (eplist)
}
\arguments{
  \item{eplist}{An \code{ep.list} object. A list of edge
  probability matrices}
}
\description{
  Plot a \code{\link[lattice]{levelplot}} displaying the
  edge probabilities.
}
\examples{
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
x <- data.frame(x1 = x1, x2 = x2)

exact <- posterior(data = x, "exact")
mcmc <- posterior(data = x, "mc3", nSamples = 500, nBurnin = 100)

myep1 <- ep(exact)
myep2 <- ep(mcmc)
if (require(lattice)){
  levelplot(ep.list(Exact = myep1, MCMC = myep2))
}
}

