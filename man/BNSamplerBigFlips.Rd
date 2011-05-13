\name{BNSamplerBigFlips}
\alias{BNSamplerBigFlips}
\title{Create a MCMC sampler for Bayesian Networks.}
\usage{BNSamplerBigFlips(data, initial, prior, return="network",
    logScoreFUN=defaultLogScoreFUN(),
    logScoreParameters=list(hyperparameters = "qi"), verbose=F,
    keepTape=F)}
\description{Create a MCMC sampler for Bayesian Networks. The sampler samples Bayesian
Networks (ie models).}
\details{.....}
\value{A function, which when called draws the next sample of the MCMC.}
\arguments{\item{data}{The data.}
\item{initial}{An object of class 'bn'. The starting value of the MCMC.}
\item{prior}{A function that returns the prior score of the supplied
bn.}
\item{return}{Either "network" or "contingency".}
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
\item{logScoreParameters}{...}
\item{verbose}{A logical of length 1, indicating whether verbose output
should be printed.}
\item{keepTape}{A logical of length 1, indicating whether a full log
('tape') of the MCMC sampler should be kept. Enabling this
option can be very memory-intensive.}}
