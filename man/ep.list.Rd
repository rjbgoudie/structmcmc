\name{ep.list}
\alias{ep.list}
\title{List of posterior edge probabilities.}
\usage{ep.list(...)}
\description{List of posterior edge probabilities.}
\details{method description}
\seealso{\code{\link{ep}}}
\arguments{\item{...}{Further arguments passed to method}}
\examples{x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
x <- data.frame(x1 = x1, x2 = x2)

exact <- posterior(data = x, "exact")
mcmc <- posterior(data = x, "mc3", nSamples = 500, nBurnin = 100)

myep1 <- ep(exact)
myep2 <- ep(mcmc)
ep.list(Exact = myep1, MCMC = myep2)}
