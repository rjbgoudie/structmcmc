\name{ep.bnpostmcmc}
\alias{ep.bnpostmcmc}
\title{Posterior edge probabiities.}
\usage{
  \method{ep}{bnpostmcmc} (x, nbin = 1, start, end,
    method = "et", verbose = F, ...)
}
\arguments{
  \item{x}{An object of class 'bnpostmcmc'}

  \item{nbin}{A numberic vector of length 1 specifying the
  number of bins into which to divide the MCMC samples. The
  edge probabilities are computed separately for each bin.}

  \item{start}{...}

  \item{end}{...}

  \item{method}{Either "et" (the default), "flatten" or
  "tabulate". If the edge totals are not available, method
  "tabulate" is used if possible. Only "flatten" is
  available if \code{nbin != 1}.}

  \item{verbose}{...}

  \item{...}{Further arguments passed to ep.parental.list()
  for method = "flatten", or ep.table() for method =
  "tabulate"}
}
\value{
  if nbin == 1: A matrix of class 'ep' with entry (i,j)
  containing the probability of an edge from node i --> j
  if nbin > 1: A list of class ep.list, containing matrices
  as described above for each of the nbin bins into which
  the parental.list was split
}
\description{
  Computes the edge probabilities implied by the MCMC
  samples contained in the 'bnpostmcmc' object x.
}
\examples{
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
dat <- data.frame(x1 = x1, x2 = x2)

prior <- function(net) 1
initial <- bn(c(), c())

sampler <- BNSampler(dat, initial, prior)
samples <- draw(sampler, n = 50)
mpost <- bnpostmcmc(sampler, samples)

ep(mpost)
}
\seealso{
  \code{\link{ep}}, \code{\link{ep.bnpostmcmc.list}}
}

