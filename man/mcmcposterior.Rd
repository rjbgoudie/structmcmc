\name{mcmcposterior}
\alias{mcmcposterior}
\title{Posterior distribution on Bayesian networks.}
\usage{mcmcposterior(data, sampler=BNSampler, prior=priorUniform(),
    logScoreFUN=logScoreMultDirFUN(),
    logScoreParameters=list(hyperparameters = "qi"), constraint,
    maxNumberParents, nSamples=50000, nBurnin=10000,
    initial=empty(ncol(data), "bn"), verbose=T)}
\description{Posterior distribution on Bayesian networks.}
\details{Use MCMC to approximate the posterior distribution}
\value{A \code{bnpostmcmc} object.}
\seealso{For more control, use the MCMC sampler directly, 
e.g. \code{\link{BNSampler}}. See also \code{\link{posterior}}.
Example priors \code{\link{priorGraph}}, \code{\link{priorUniform}}.}
\arguments{\item{data}{The data.}
\item{sampler}{A BNSampler. eg BNSampler or BNGibbsSampler etc}
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
