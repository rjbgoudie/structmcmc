\name{logScoreMultDirFUN}
\alias{logScoreMultDirFUN}
\title{Multinomial-Dirichlet Log marginal likelihood.}
\usage{
  logScoreMultDirFUN()
}
\value{
  A list consisting of the functions that perform each of
  the following roles. \describe{ \item{offline}{A function
  that computes the logScore of a Bayesian Network}
  \item{online}{A function that incrementally computes the
  logScore of a Bayesian Network} \item{local}{A function
  that computes the local logScore of a Bayesian Network}
  \item{prepare}{A function that prepares the data, and any
  further pre-computation required by the logScore
  functions.} }
}
\description{
  This function returns a list of functions that are used
  for MCMC computation for Multinomial-Dirichlet models.
}
\seealso{
  \code{\link{logScoreNormalFUN}},
  \code{\link{logScoreMultDir}},
  \code{\link{logScoreNormal}}
}

