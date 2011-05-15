\name{map.bnpostmcmc}
\alias{map.bnpostmcmc}
\title{Maximum aposteriori graph.}
\usage{\method{map}{bnpostmcmc} (x, ...)}
\description{Maximum aposteriori graph.}
\details{Returns the most commonly encountered graph(s) during the MCMC sampling
'x'. ie the maximum aposteriori graph(s).}
\value{EITHER an object of class 'bn' containing the most commonly
encountered graph,
OR an object of class 'bn.list' containing a list of the equally most-
commonly encountered graphs (if two or more were equally most-
commonly) encountered}
\seealso{\code{\link{map}}}
\arguments{\item{x}{An object of class 'bnpostmcmc'}
\item{...}{Further arguments (unused)}}
