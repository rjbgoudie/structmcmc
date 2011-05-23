\name{xyplot.cumtvd}
\alias{xyplot.cumtvd}
\title{Plot cumulative total variation distance.}
\usage{\method{xyplot}{cumtvd} (cumtvd)}

\description{
  Plot cumulative total variation distance.
}

\details{
  method description
}
\seealso{\code{\link{cumtvd}}}
\arguments{
  \item{cumtvd}{...}
}
\examples{x1 <- factor(c(1, 1, 0, 1, 1, 1))
x2 <- factor(c(0, 1, 0, 1, 0, 1))
x3 <- factor(c(0, 1, 0, 0, 0, 0))
x <- data.frame(x1 = x1, x2 = x2, x3 = x3)

exact <- posterior(x, "exact")
exactgp <- gp(exact)
mcmc1 <- posterior(x, "mc3", nSamples = 1000, nBurnin = 100)
mcmc2 <- posterior(x, "mc3", nSamples = 1000, nBurnin = 100)

tvd <- cumtvd(exactgp = exactgp,
bnpostmcmclist = bnpostmcmc.list(mcmc1, mcmc2))
if (require(lattice)){
xyplot(tvd)
}}
