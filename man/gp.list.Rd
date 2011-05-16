\name{gp.list}
\alias{gp.list}
\title{List of posterior graph probabilities.}
\usage{gp.list(...)}
\description{List of posterior graph probabilities.}
\details{method description}
\arguments{\item{...}{Further arguments passed to method}}
\examples{x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
x <- data.frame(x1 = x1, x2 = x2)

exact <- posterior(data = x, "exact")
mcmc <- posterior(data = x, "mc3", nSamples = 500, nBurnin = 100)

mygp1 <- gp(exact)
mygp2 <- gp(mcmc)
gp.list(Exact = mygp1, MCMC = mygp2)}
