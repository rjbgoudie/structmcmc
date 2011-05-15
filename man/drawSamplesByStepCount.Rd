\name{drawSamplesByStepCount}
\alias{drawSamplesByStepCount}
\title{Draw samples from a MCMC sampler, by step count.}
\usage{drawSamplesByStepCount(sampler, n, burnin, burninleft, lengthout, thin,
    verbose)}
\description{Draw samples from a MCMC sampler, by step count.}
\details{Draws a specific number of samples from an MCMC sampler.}
\seealso{\code{\link{draw}}, \code{\link{drawSamplesByTime}}}
\arguments{\item{sampler}{A sampler object}
\item{n}{The number of samples to draw.}
\item{burnin}{The number of samples to discard from the beginning of
the sample (includes previous runs)}
\item{burninleft}{The number of samples still to be discarded.}
\item{lengthout}{The length of the final samples object}
\item{thin}{The frequency with which samples should be kept. eg for
\code{thin = 3}, every third sample will be kept.}
\item{verbose}{A logical. Should a progress bar be displayed?}}
\examples{set.seed(310)
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
dat <- data.frame(x1 = x1, x2 = x2)

prior <- function(net) 1
initial <- bn(c(), c())

sampler <- BNSampler(dat, initial, prior)
samples <- draw(sampler, n = 100)}
