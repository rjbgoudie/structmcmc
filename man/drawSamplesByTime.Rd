\name{drawSamplesByTime}
\alias{drawSamplesByTime}
\title{Draw samples from a MCMC sampler, by time.}
\usage{
  drawSamplesByTime(sampler, time, burnin, burninleft,
    thin, verbose, samplesIncrement = 10000)
}
\arguments{
  \item{sampler}{A sampler object}

  \item{time}{The number of seconds to spend drawing
  samples.}

  \item{burnin}{The number of samples to discard from the
  beginning of the sample (includes previous runs)}

  \item{burninleft}{The number of samples still to be
  discarded.}

  \item{thin}{The frequency with which samples should be
  kept. eg for \code{thin = 3}, every third sample will be
  kept.}

  \item{verbose}{A logical. Should a progress bar be
  displayed?}

  \item{samplesIncrement}{The size by which the samples
  should be incremented.}
}
\description{
  Draws samples from an MCMC sampler for a specified number
  of seconds.
}
\examples{
set.seed(310)
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
dat <- data.frame(x1 = x1, x2 = x2)

prior <- function(net) 1
initial <- bn(c(), c())

sampler <- BNSampler(dat, initial, prior)
samples <- draw(sampler, time = 15)
}
\seealso{
  \code{\link{draw}}, \code{\link{drawSamplesByStepCount}}
}

