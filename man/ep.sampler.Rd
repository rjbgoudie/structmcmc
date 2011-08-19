\name{ep.sampler}
\alias{ep.sampler}
\title{Extract posterior edge probabiities from a sampler.}
\usage{
  \method{ep}{sampler} (x, start, end, verbose = F, ...)
}
\arguments{
  \item{x}{An MCMC sampler 'function'}

  \item{start}{...}

  \item{end}{...}

  \item{verbose}{...}

  \item{...}{Further arguments passed to ep.parental.list()
  for method = "flatten", or ep.table() for method =
  "tabulate"}
}
\value{
  A matrix of class 'ep' with entry (i,j) containing the
  probability of an edge from node i --> j
}
\description{
  Computes the edge probabilities implied by the MCMC
  samples contained in a MCMC sampler function.
}
\examples{
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
dat <- data.frame(x1 = x1, x2 = x2)

prior <- function(net) 1
initial <- bn(c(), c())

sampler <- BNSampler(dat, initial, prior)
samples <- draw(sampler, n = 50)

ep(sampler)
}
\seealso{
  \code{\link{ep}}, \code{\link{ep.bnpostmcmc.list}}
}

