\name{statistics.sampler}
\alias{statistics.sampler}
\title{Extract statistics from a sampler.}
\usage{
  statistics.sampler(x, names, ...)
}
\arguments{
  \item{x}{A sampler}

  \item{names}{Which statistics to extract? A character
  vector of statistics collected during the MCMC run.}

  \item{...}{Further arguments passed to method}
}
\value{
  A list of statistics collected during the MCMC run.
}
\description{
  Extracts the statistics collected during an MCMC run
}
\seealso{
  \code{\link{statistics}}
}

