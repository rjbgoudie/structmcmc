\name{posterior}
\alias{posterior}
\title{Posterior distribution on Bayesian networks...}
\usage{posterior(data, method="mh-mcmc", prior=function(x) 1,
    logScoreFUN=defaultLogScoreFUN(),
    logScoreParameters=list(hyperparameters = "qi"), constraint,
    maxNumberParents=Inf, nSamples=50000, nBurnin=10000, initial,
    verbose=F)}
\description{Posterior distribution on Bayesian networks}
\details{Use one of a number of methods to get the posterior distribution}
\value{Either a \code{bnpost} or a \code{bnpostmcmc} object.}
\arguments{\item{data}{The data.}
\item{method}{One of "exact", "mh-mcmc", "gibbs", "mj-mcmc"}
\item{prior}{A function that returns the prior score of the
supplied bn.}
\item{logScoreFUN}{A list of three elements:
offline: A function that computes the logScore
of a Bayesian Network
online:  A function that incrementally computes
the logScore of a Bayesian Network
prepare: A function that prepares the data, and
any further pre-computation required by
the logScore functions.}
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
