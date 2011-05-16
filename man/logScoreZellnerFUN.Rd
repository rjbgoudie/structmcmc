\name{logScoreZellnerFUN}
\alias{logScoreZellnerFUN}
\title{Normal Log marginal likelihood.}
\description{Normal Log marginal likelihood.}
\details{This function returns a list of functions that are used for MCMC
computation for Normal models, with Zellner g-priors.}
\value{A list consisting of the functions that perform each of the 
following roles.
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
\seealso{\code{\link{logScoreMultDirFUN}}, \code{\link{logScoreMultDir}},
\code{\link{logScoreZellner}}}
