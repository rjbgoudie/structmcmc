\name{bnpost}
\alias{bnpost}
\title{BN Exact Posterior.}
\usage{
  bnpost(bnspace, logScore, data, logScoreFUN =
  logScoreMultDir)
}
\arguments{
  \item{bnspace}{...}

  \item{logScore}{...}

  \item{data}{...}

  \item{logScoreFUN}{A list of four elements: \describe{
  \item{offline}{A function that computes the logScore of a
  Bayesian Network} \item{online}{A function that
  incrementally computes the logScore of a Bayesian
  Network} \item{local}{A function that computes the local
  logScore of a Bayesian Network} \item{prepare}{A function
  that prepares the data, and any further pre-computation
  required by the logScore functions.} } For
  Multinomial-Dirichlet models,
  \code{\link{logScoreMultDirFUN}} returns the appropriate
  list; for Normal models with Zellner g-priors,
  \code{\link{logScoreNormalFUN}} returns the appropriate
  list.}
}
\description{
  method description
}
\examples{
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
dat <- data.frame(x1 = x1, x2 = x2)

bnspace <- enumerateBNSpace(2)
lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)

gp(post)
ep(post)
}
\seealso{
  \code{\link{bnpostmcmc}}, \code{\link{top.bnpost}},
  \code{\link{map.bnpost}}, \code{\link{gp.bnpost}},
  \code{\link{entropy.bnpost}}, \code{\link{tp.bnpost}},
  \code{\link{eht.bnpost}}
}

