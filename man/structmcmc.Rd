\name{structmcmc}
\title{The structmcmc package.}
\description{The structmcmc package.}
\docType{package}
\alias{structmcmc}
\alias{package-structmcmc}
\seealso{\code{\link{posterior}}}
\examples{# Setup data frame
x1 <- factor(c("a", "a", "g", "c", "c", "a", "g", "a", "a"))
x2 <- factor(c(2, 2, 4, 3, 1, 4, 4, 4, 1))
x3 <- factor(c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
x <- data.frame(x1 = x1, x2 = x2, x3 = x3)

# Draw samples from the posterior using MC^3.
set.seed(1234)
mcmc <- posterior(data = x, method = "mh-mcmc",
nSamples = 1000, nBurnin = 100)

# Compute and plot estimated edge probabilities
epmcmc <- ep(mcmc)

# Exact evaluation by exhaustive enumeration
exact <- posterior(x, "exact")
epexact <- ep(exact)

# Comparing multiple MCMC runs
mcmc2 <- posterior(x, "mh-mcmc", nSamples = 1000, nBurnin = 100)
epmcmc2 <- ep(mcmc2)}
