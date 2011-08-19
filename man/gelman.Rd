\name{gelman}
\alias{gelman}
\title{Gelman and Rubin's convergence diagnostic}
\usage{
  gelman(x, ...)
}
\arguments{
  \item{x}{...}

  \item{...}{Further arguments passed to method}
}
\description{
  The 'potential scale reduction factor' is calculated for
  each variable in \code{x}, together with upper and lower
  confidence limits. Computation is performed using
  \code{\link[coda]{gelman.diag}}.
}
\seealso{
  \code{\link{gelman.samplers}}. Computation performed
  using \code{\link[coda]{gelman.diag}}.
}

