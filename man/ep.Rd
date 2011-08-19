\name{ep}
\alias{ep}
\title{Posterior edge probabilities.}
\usage{
  ep(x, ...)
}
\arguments{
  \item{x}{An object}

  \item{...}{Further arguments, passed to method}
}
\description{
  Compute the posterior edge probabilities, as given by the
  supplied object.
}
\examples{
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
dat <- data.frame(x1 = x1, x2 = x2)

initial <- bn(c(), c())
prior <- priorUniform(initial)

sampler <- BNSampler(dat, initial, prior)
samples <- draw(sampler, n = 50)
mpost <- bnpostmcmc(sampler, samples)

ep(mpost)

initial <- bn(c(), c(1))
sampler2 <- BNSampler(dat, initial, prior)
samples2 <- draw(sampler2, n = 50)
mpost2 <- bnpostmcmc(sampler2, samples2)

ep(bnpostmcmc.list(mpost, mpost2))
}
\seealso{
  \code{\link{ep.bnpostmcmc.list}},
  \code{\link{ep.parental.list}}, \code{\link{ep.bnpost}},
  \code{\link{ep.table}},
  \code{\link{ep.parental.contingency}}
}

