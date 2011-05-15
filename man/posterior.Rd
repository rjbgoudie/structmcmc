\name{posterior}
\alias{posterior}
\title{Posterior distribution on Bayesian networks.}
\usage{posterior(data, method="mc3", prior=function(x) 1,
    logScoreFUN=defaultLogScoreFUN(),
    logScoreParameters=list(hyperparameters = "qi"), constraint,
    maxNumberParents=Inf, nSamples=50000, nBurnin=10000, initial,
    verbose=F)}
\description{Posterior distribution on Bayesian networks.}
\details{Use one of a number of methods to get the posterior distribution.}
\value{Either a \code{bnpost} or a \code{bnpostmcmc} object.}
\seealso{For more control, use the MCMC sampler directly, 
e.g. \code{\link{BNSampler}}.}
\arguments{\item{data}{The data.}
\item{method}{One of "exact", "mc3", "gibbs", "mj-mcmc". "mh-mcmc" is a 
synonym of "mc3".}
\item{prior}{A function that returns the prior score of the
supplied bn.}
\item{logScoreFUN}{A list of four elements:
\describe{
\item{offline}{A function that computes the logScore of a Bayesian 
Network}
\item{online}{A function that incrementally computes the logScore of a 
Bayesian Network}
\item{local}{A function that computes the local logScore of a 
Bayesian Network}
\item{prepare}{A function that prepares the data, and any further 
pre-computation required by the logScore functions.}
}}
\item{logScoreParameters}{A list of parameters that are passed to
logScoreFUN.}
\item{constraint}{A matrix of dimension ncol(data) x ncol(data) giving
constraints to the sample space.
The (i, j) element is
1  if the edge i -> j is required
-1 if the edge i -> is excluded.
0  if the edge i -> j is not constrained.
The diagonal of constraint must be all 0.}
\item{maxNumberParents}{Integer of length 1. The maximum number of
parents of any node.}
\item{nSamples}{The number of samples to be draw (only applies to MCMC)}
\item{nBurnin}{The number of samples to discard from the beginning of
the sample.}
\item{initial}{An object of class 'bn'. The starting value of the
MCMC.}
\item{verbose}{A logical. Should a progress bar be displayed?}}
\examples{x1 <- factor(c("a", "a", "g", "c", "c", "a", "g", "a", "a"))
x2 <- factor(c(2, 2, 4, 3, 1, 4, 4, 4, 1))
x3 <- factor(c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
x <- data.frame(x1 = x1, x2 = x2, x3 = x3)

set.seed(1234)
mcmc <- posterior(data = x, "mc3", nSamples = 500, nBurnin = 100)
ep(mcmc)}
