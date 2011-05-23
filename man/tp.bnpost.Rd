\name{tp.bnpost}
\alias{tp.bnpost}
\title{Transition probabilities.}
\usage{\method{tp}{bnpost} (x, sampler="mh", allowFlips=T,
    verbose=F, ...)}

\description{
  Transition probabilities.
}

\details{
  Computes the matrix transition probabilities for the
  specified sampler.
}
\value{A matrix of transition probabilities.}
\seealso{\code{\link{tp}}}
\arguments{
  \item{x}{An object of class "bnpost"}
  \item{sampler}{Which sampler to use. Only "mh" for Metropolis-Hastings
is implemented}
  \item{allowFlips}{A logical of length 1, specifying whether the sampler is
allowed to reverse the direction of single edges?}
  \item{verbose}{...}
  \item{...}{Further arguments (unused)}
}
\examples{x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
dat <- data.frame(x1 = x1, x2 = x2)

bnspace <- enumerateBNSpace(2)
lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)

tp(post)}
