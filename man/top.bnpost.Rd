\name{top.bnpost}
\alias{top.bnpost}
\title{Top graph from BN Posterior.}
\usage{
  \method{top}{bnpost} (x, head = 10, ...)
}
\arguments{
  \item{x}{...}

  \item{head}{...}

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
lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "bdeu")
post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)

top(post)
}
\seealso{
  \code{\link{top}}, \code{\link{map.bnpost}}
}

