\name{top.bnpostmcmc}
\alias{top.bnpostmcmc}
\title{Top graph from BN Posterior.}
\usage{\method{top}{bnpostmcmc} (x, head=10, ...)}
\description{Top graph from BN Posterior.}
\details{Returns the most commonly encountered graphs during the MCMC sampling
'x'. The top 'head' graphs with respect to MCMC sampling are returned.
ie if the MCMC sampler has converged, the top graphs with respect to the
posterior distribution on graphs will be returned.}
\value{if head == 1:
EITHER an object of class 'bn' containing the most commonly
encountered graph,
OR an object of class 'bn.list' containing a list of the equally most-
commonly encountered graphs (if two or more were equally most-
commonly) encountered
if head > 1:
an object of class 'bn.list' containing a list of the equally most-
commonly encountered graphs}
\arguments{\item{x}{An object of class 'bnpostmcmc'}
\item{head}{The top head graphs will be returned. If the 'head'th
most-commonly encountered graph ties in frequency with other
graphs, all of the ties will be returned.}
\item{...}{Further arguments (unused)}}
