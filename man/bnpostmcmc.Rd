\name{bnpostmcmc}
\alias{bnpostmcmc}
\title{BN Posterior from MCMC.}
\usage{bnpostmcmc(sampler, samples, logScoreFUN)}
\description{BN Posterior from MCMC.}
\details{method description}
\arguments{\item{sampler}{...}
\item{samples}{...}
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
}}}
