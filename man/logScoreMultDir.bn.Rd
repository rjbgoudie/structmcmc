\name{logScoreMultDir.bn}
\alias{logScoreMultDir.bn}
\title{Compute the log marginal likelihood of the supplied Bayesian Network.}
\usage{\method{logScoreMultDir}{bn} (x, data, cache=new.env(hash = T),
    hyperparameters="qi", checkInput=T, ...)}
\description{Compute the log marginal likelihood of the supplied Bayesian Network.}
\details{The data must be discrete. The conditional distributions of each
random variable, conditional on its parents are assumed to be
multinomial, with Dirichlet priors for the parameters.

The notation here roughly follows
Mukherjee and Speed (2008) Network inference using informative priors.
PNAS 105 (38) 14313-14318, doi: 10.1073/pnas.0802272105}
\value{A numeric vector of length 1, giving the log marginal likelihood.
The environment 'cache' will also be updated because its scope is
global.}
\seealso{\code{\link{logScoreMultDir}},
\code{\link{logScoreMultDir.bn.list}}}
\arguments{\item{x}{An object of class "bn". The Bayesian Network by for which the
marginal likelihood is computed.}
\item{data}{A data.frame, with columns being factors giving the values of
each random variable.}
\item{cache}{Optionally, provide an environment with cached local scores
for this data.}
\item{hyperparameters}{A character vector of length one. Either "qi",
"one", or "point9"}
\item{checkInput}{A logical of length 1, specifying whether to check the
inputs to the function.}
\item{...}{Further arguments, currently unused}}
