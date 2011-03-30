\name{logScoreMultDirOffline}
\alias{logScoreMultDirOffline}
\title{Compute the log marginal likelihood of the supplied Bayesian Network.}
\usage{logScoreMultDirOffline(x, logScoreParameters, cache=new.env(hash = T),
    checkInput=T)}
\description{Compute the log marginal likelihood of the supplied Bayesian Network.}
\details{This function is an alternative interface to logScoreMultDir.
This interface is required by the MCMC sampler.}
\value{A numeric vector of length 1, giving the log marginal likelihood.
The environment 'cache' will also be updated because its scope is
global.}
\arguments{\item{x}{An object of class "bn". The Bayesian Network by for
which the marginal likelihood is computed.}
\item{logScoreParameters}{...}
\item{cache}{Optionally, provide an environment with cached
local scores for this data.}
\item{checkInput}{A logical of length 1, specifying whether to check
the inputs to the function.}}
