\name{eht.bnpost}
\alias{eht.bnpost}
\title{Expected hitting time.}
\usage{
  \method{eht}{bnpost} (x, head = 5, tp = NULL, ...)
}
\arguments{
  \item{x}{An object of class "bnpost".}

  \item{head}{The expected hitting time from the top 'head'
  ranked graphs, with respect to the posterior will be
  returned. All are computed, so this setting has no effect
  on the run-time.}

  \item{tp}{Optionally provide the transition probability
  matrix.}

  \item{...}{Further arguments (unused)}
}
\value{
  A vector of the expected hitting times.
}
\description{
  Computes the expected hitting times to the posterior
  modal graph from the top 'head' other graphs.
}
\details{
  See Norris (1998) Markov Chains. Cambridge. Theorem
  1.3.5, p17.
}
\seealso{
  \code{\link{eht}}, \code{\link{tp}}, \code{\link{hp}}
}

