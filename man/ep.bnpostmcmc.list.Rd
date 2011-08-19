\name{ep.bnpostmcmc.list}
\alias{ep.bnpostmcmc.list}
\title{Posterior edge probabilities.}
\usage{
  \method{ep}{bnpostmcmc.list} (x, start, end, ...)
}
\arguments{
  \item{x}{...}

  \item{start}{...}

  \item{end}{...}

  \item{...}{further arguments}
}
\description{
  method description
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

initial <- bn(c(), c(1))
sampler2 <- BNSampler(dat, initial, prior)
samples2 <- draw(sampler2, n = 50)
mpost2 <- bnpostmcmc(sampler2, samples2)

ep(bnpostmcmc.list(mpost, mpost2))
}
\seealso{
  \code{\link{ep}}, \code{\link{ep.bnpostmcmc}}
}

