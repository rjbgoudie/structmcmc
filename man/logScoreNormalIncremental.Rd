\name{logScoreNormalIncremental}
\alias{logScoreNormalIncremental}
\title{Normal-inverse-gamma (with g-prior) Log marginal likelihood (online).}
\usage{logScoreNormalIncremental(currentBN, proposalBN, heads,
    logScoreParameters, cache, checkInput=T)}

\description{
  Normal-inverse-gamma (with g-prior) Log marginal
  likelihood (online).
}

\details{
  Compute the difference in log marginal likelihood of the
  supplied Bayesian Networks, quickly.

  This is a fast, incremental version of logScoreNormal.

  The data is scored as continuous, using a form of the
  Zellner Prior.
}
\value{logscore(proposalBN) - logscore(currentBN)}
\seealso{\code{\link{logScoreNormal}},
\code{\link{logScoreNormalOffline}}}
\arguments{
  \item{currentBN}{An object of class "bn".}
  \item{proposalBN}{An object of class "bn".}
  \item{heads}{A numeric vector, specifying which nodes have different
parents in currentBN and proposalBN.}
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
