\name{draw}
\alias{draw}
\title{Draw samples from a MCMC sampler.}
\usage{
  draw(sampler, n = F, time = F, burnin = 0, thin = 1,
  verbose = T)
}
\arguments{
  \item{sampler}{A sampler object}

  \item{n}{The number of samples to draw. Set this to
  \code{FALSE} if using the \code{time} argument.}

  \item{time}{The number of seconds to spend drawing
  samples. Set this to \code{FALSE} if using the \code{n}
  argument.}

  \item{burnin}{The number of samples to discard from the
  beginning of the sample.}

  \item{thin}{The frequency with which samples should be
  kept. eg for \code{thin = 3}, every third sample will be
  kept.}

  \item{verbose}{A logical. Should a progress bar be
  displayed?}
}
\description{
  Draws samples from an MCMC sampler. The length of the run
  can be specified either by the number of samples to be
  drawn, or the length of time that the sampler runs.
}
\examples{
set.seed(310)
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
dat <- data.frame(x1 = x1, x2 = x2)

prior <- function(net) 1
initial <- bn(c(), c())

sampler <- BNSampler(dat, initial, prior)
samples <- draw(sampler, n = 5)
}
\seealso{
  \code{\link{drawSamplesByTime}},
  \code{\link{drawSamplesByStepCount}}.
  \code{\link{BNSampler}}, \code{\link{BNGibbsSampler}},
  \code{\link{BNSamplerMJ}}
}

