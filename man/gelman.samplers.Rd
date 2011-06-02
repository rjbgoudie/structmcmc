\name{gelman.samplers}
\alias{gelman.samplers}
\title{Gelman and Rubin's convergence diagnostic...}
\usage{\method{gelman}{samplers} (x, names="nEdges", transform,
    ...)}

\description{
  Gelman and Rubin's convergence diagnostic
}

\details{
  The 'potential scale reduction factor' is calculated for
  each variable in \code{x}, together with upper and lower
  confidence limits. Computation is performed using
  \code{\link[coda]{gelman.diag}}.
}
\value{An object of class \code{gelman.diag}}
\seealso{Computation performed using \code{\link[coda]{gelman.diag}}.}
\arguments{
  \item{x}{An MCMC sampler}
  \item{names}{A character vector of statistics collected during the
MCMC run, for which the 'potential scale reduction factor' should be
computed.}
  \item{transform}{A list of functions (or a single function), which will
be applied to the statistics. This can be used to improve the
normality of the statistics, which is a requirement for Gelman-Rubin's
statistic to be meaningful. The default \code{NULL} value applies no
transform to any statistic (equivalent to \code{\link{identity}})}
  \item{...}{Further arguments, currently unused}
}
