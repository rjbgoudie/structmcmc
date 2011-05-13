\name{exactposterior}
\alias{exactposterior}
\title{Posterior distribution on Bayesian networks...}
\usage{exactposterior(data, prior, logScoreFUN, logScoreParameters,
    constraint, maxNumberParents, verbose)}
\description{Posterior distribution on Bayesian networks}
\details{Use one of a number of methods to get the posterior distribution}
\value{A \code{bnpost} object.}
\arguments{\item{data}{The data.}
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
\item{verbose}{A logical. Should a progress bar be displayed?}}
