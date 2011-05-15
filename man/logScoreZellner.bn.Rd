\name{logScoreZellner.bn}
\alias{logScoreZellner.bn}
\title{Normal Log marginal likelihood.}
\usage{\method{logScoreZellner}{bn} (x, data, cache=new.env(hash = T),
    checkInput=T, ...)}
\description{Normal Log marginal likelihood.}
\details{Compute the log marginal likelihood of the supplied Bayesian Network.

The data is scored as continuous, using a form of the Zellner Prior.}
\value{A numeric vector of length 1, giving the log marginal likelihood.
The environment 'cache' will also be updated because its scope is
global.}
\arguments{\item{x}{An object of class "bn". The Bayesian Network by
for which the marginal likelihood is computed.}
\item{data}{A matrix with columns giving the value of each random
variable.}
\item{cache}{Optionally, provide an environment with cached
local scores for this data.}
\item{checkInput}{A logical of length 1, specifying whether to check
the inputs to the function.}
\item{...}{Further arguments (unused)}}
