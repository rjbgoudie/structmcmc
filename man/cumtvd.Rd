\name{cumtvd}
\alias{cumtvd}
\title{Cumulative total variation distance.}
\usage{
  cumtvd(exactgp, bnpostmcmclist, start = 1, end,
    nbin = (end - start + 1)/100)
}
\arguments{
  \item{exactgp}{...}

  \item{bnpostmcmclist}{...}

  \item{start}{...}

  \item{end}{...}

  \item{nbin}{...}

  \item{...}{Further arguments passed to method}
}
\description{
  method description
}
\examples{
x1 <- factor(c(1, 1, 0, 1, 1, 1))
x2 <- factor(c(0, 1, 0, 1, 0, 1))
x3 <- factor(c(0, 1, 0, 0, 0, 0))
x <- data.frame(x1 = x1, x2 = x2, x3 = x3)

exact <- posterior(x, "exact")
exactgp <- gp(exact)
mcmc1 <- posterior(x, "mc3", nSamples = 1000, nBurnin = 100)
mcmc2 <- posterior(x, "mc3", nSamples = 1000, nBurnin = 100)

tvd <- cumtvd(exactgp = exactgp,
              bnpostmcmclist = bnpostmcmc.list(mcmc1, mcmc2))
}
\seealso{
  \code{\link{xyplot.cumtvd}}, \code{\link{cumtvd.list}}
}

