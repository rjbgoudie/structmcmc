\name{localLogScoreMultDir}
\alias{localLogScoreMultDir}
\title{Local Multinomial-Dirichlet Log marginal likelihood.}
\usage{
  localLogScoreMultDir(node, parents, logScoreParameters,
    cache, checkInput = T)
}
\arguments{
  \item{node}{A numeric vector of length 1. The node to
  compute the local log score for.}

  \item{parents}{A numeric vector. The parents of node.}

  \item{logScoreParameters}{A list with the following
  components: \describe{ \item{data}{A matrix (NOT
  data.frame), with columns being integers in the range 0,
  1, 2, ....  giving the values of each random variable.
  **** The integers MUST start numbering at 0 NOT 1 ****}
  \item{nl}{A numeric vector of length ncol(data),
  specifying the number of levels that each random variable
  takes.} \item{hyperparameters}{A character vector of
  length one.  Either "qi", "one", or "point9"} }}

  \item{cache}{Optionally, provide an environment with
  cached local scores for this data.}

  \item{checkInput}{A logical of length 1, specifying
  whether to check the inputs to the function.}
}
\value{
  A numeric vector of length 1, giving the log marginal
  likelihood.  The environment 'cache' will also be updated
  because its scope is global.
}
\description{
  Compute the LOCAL log marginal likelihood of the supplied
  Bayesian Networks. ie the contribution to the log
  marginal liklihood from one individual node.
}
\details{
  The data must be discrete. The conditional distributions
  of each random variable, conditional on its parents are
  assumed to be multinomial, with Dirichlet priors for the
  parameters.

  The notation here roughly follows Mukherjee and Speed
  (2008) Network inference using informative priors. PNAS
  105 (38) 14313-14318, doi: 10.1073/pnas.0802272105
}
\seealso{
  \code{\link{logScoreMultDir}},
  \code{\link{logScoreMultDirIncremental}},
  \code{\link{logScoreMultDirOffline}}
}

