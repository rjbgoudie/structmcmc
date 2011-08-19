\name{plot.bnpostmcmc}
\alias{plot.bnpostmcmc}
\title{Plot top graphs.}
\usage{
  \method{plot}{bnpostmcmc} (x, top = NULL, ...)
}
\arguments{
  \item{x}{An object of class "bnpostmcmc" or "bnpost"}

  \item{top}{Optionally provide pre-computed top(x)}

  \item{...}{Further arguments, passed to
  \code{\link[parental]{grplot}}}
}
\value{
  A plot of the top graph, with their marginal likelihoods
  (without priors)
}
\description{
  Plot the the top() Bayesian Networks from a posterior
  distribution. The top graphs are those with the highest
  score with respect to the posterior distribution, which
  for converged MCMC will be most commonly encountered
  graph.
}
\seealso{
  \code{\link{levelplot.bnpostmcmc}},
  \code{\link{bnpostmcmc}}, \code{\link{bnpost}}
}

