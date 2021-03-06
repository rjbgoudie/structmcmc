\name{logScoreMultDir.bn}
\alias{logScoreMultDir.bn}
\title{Compute the log marginal likelihood of the supplied Bayesian Network.}
\usage{
  \method{logScoreMultDir}{bn} (x, data,
    cache = new.env(hash = T), hyperparameters = "bdeu",
    checkInput = T, ...)
}
\arguments{
  \item{x}{An object of class "bn". The Bayesian Network by
  for which the marginal likelihood is computed.}

  \item{data}{A data.frame, with columns being factors
  giving the values of each random variable.}

  \item{cache}{Optionally, provide an environment with
  cached local scores for this data.}

  \item{hyperparameters}{A character vector of length one.
  Either "bdeu", "qi", "one", or "point9"}

  \item{checkInput}{A logical of length 1, specifying
  whether to check the inputs to the function.}

  \item{...}{Further arguments, currently unused}
}
\value{
  A numeric vector of length 1, giving the log marginal
  likelihood.  The environment 'cache' will also be updated
  because its scope is global.
}
\description{
  The data must be discrete. The conditional distributions
  of each random variable, conditional on its parents are
  assumed to be multinomial, with Dirichlet priors for the
  parameters.
}
\details{
  The notation here roughly follows Mukherjee and Speed
  (2008) Network inference using informative priors. PNAS
  105 (38) 14313-14318, doi: 10.1073/pnas.0802272105
}
\examples{
x1 <- factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
x2 <- factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
x3 <- factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
data <- data.frame(x1 = x1, x2 = x2,  x3 = x3)
logScoreMultDir(bn(c(), c(1), c(2)), data)
}
\seealso{
  \code{\link{logScoreMultDir}},
  \code{\link{logScoreMultDir.bn.list}}
}

