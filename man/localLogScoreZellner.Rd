\name{localLogScoreZellner}
\alias{localLogScoreZellner}
\title{Local Normal Log marginal likelihood.}
\usage{localLogScoreZellner(node, parents, logScoreParameters,
    cache, checkInput=T)}

\description{
  Local Normal Log marginal likelihood.
}

\details{
  Compute the LOCAL log marginal likelihood of the supplied
  Bayesian Networks. ie the contribution to the log
  marginal liklihood from one individual node.

  The data is scored as continuous, using a form of the
  Zellner Prior.
}
\value{A numeric vector of length 1, giving the log marginal likelihood.
The environment 'cache' will also be updated because its scope is
global.}
\seealso{\code{\link{logScoreZellner}},
\code{\link{logScoreZellnerOffline}},
\code{\link{logScoreZellnerIncremental}}.}
\arguments{
  \item{node}{A numeric vector of length 1. The node to compute the local
log score for.}
  \item{parents}{A numeric vector. The parents of node.}
  \item{logScoreParameters}{A list with the following components:
\describe{
\item{data}{A matrix with columns giving the values of each random
variable.}
\item{nl}{A numeric vector of length nNodes(currentBN), specifying the
number of levels that each random variable takes.}
}}
  \item{cache}{Optionally, provide an environment with cached local scores
for this data.}
  \item{checkInput}{A logical of length 1, specifying whether to check the
inputs to the function.}
}
