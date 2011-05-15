\name{ep.bnpostmcmc}
\alias{ep.bnpostmcmc}
\title{Posterior edge probabiities.}
\usage{\method{ep}{bnpostmcmc} (x, nbin=1, start, end, method="flatten",
    verbose=F, ...)}
\description{Posterior edge probabiities.}
\details{Computes the edge probabilities implied by the MCMC samples
contained in the 'bnpostmcmc' object x.}
\value{if nbin == 1:
A matrix of class 'ep' with entry (i,j) containing the probability of
an edge from node i --> j
if nbin > 1:
A list of class ep.list, containing matrices as described above for
each of the nbin bins into which the parental.list was split}
\arguments{\item{x}{An object of class 'bnpostmcmc'}
\item{nbin}{A numberic vector of length 1 specifying the number of bins
into which to divide the MCMC samples. The edge probabilities
are computed separately for each bin.}
\item{start}{...}
\item{end}{...}
\item{method}{Either "flatten" (the default) or "tabulate".}
\item{verbose}{...}
\item{...}{Further arguments passed to ep.parental.list() for method =
"flatten", or ep.table() for method = "tabulate"}}
