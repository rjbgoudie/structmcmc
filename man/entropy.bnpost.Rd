\name{entropy.bnpost}
\alias{entropy.bnpost}
\title{Entropy.}
\usage{
  \method{entropy}{bnpost} (x, logNetworkPriors, ...)
}
\arguments{
  \item{x}{...}

  \item{logNetworkPriors}{...}

  \item{...}{Further arguments (unused)}
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

entropy(post)
}
\seealso{
  \code{\link{entropy}}
}

