\name{logScoreMultDirOffline}
\alias{logScoreMultDirOffline}
\title{Multinomial-Dirichlet Log marginal likelihood (offline).}
\usage{logScoreMultDirOffline(x, logScoreParameters, cache=new.env(hash = T),
    checkInput=T)}
\description{Multinomial-Dirichlet Log marginal likelihood (offline).}
\details{Compute the log marginal likelihood of the supplied Bayesian Network.

This function is an alternative interface to logScoreMultDir.
This interface is required by the MCMC sampler.}
\value{A numeric vector of length 1, giving the log marginal likelihood.
The environment 'cache' will also be updated because its scope is
global.}
\seealso{\code{\link{logScoreMultDir}},
\code{\link{logScoreMultDirIncremental}}}
\arguments{\item{x}{An object of class "bn". The Bayesian Network by for which the 
marginal likelihood is computed.}
\item{logScoreParameters}{A list with the following components:
\describe{
\item{data}{A matrix (NOT data.frame), with columns being integers
in the range 0, 1, 2, ....  giving the values of each
random variable.
**** The integers MUST start  numbering at 0 NOT 1    ****}
\item{nl}{A numeric vector of length ncol(data), specifying the number 
of levels that each random variable takes.}
\item{hyperparameters}{A character vector of length one.
Either "qi", "one", or "point9"}
}}
\item{cache}{Optionally, provide an environment with cached local scores 
for this data.}
\item{checkInput}{A logical of length 1, specifying whether to check the
inputs to the function.}}
