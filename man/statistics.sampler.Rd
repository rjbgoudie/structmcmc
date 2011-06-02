\name{statistics.sampler}
\alias{statistics.sampler}
\title{Extract statistics from a sampler.}
\usage{statistics.sampler(x, names, ...)}

\description{
  Extract statistics from a sampler.
}

\details{
  Extracts the statistics collected during an MCMC run
}
\value{A list of statistics collected during the MCMC run.}
\seealso{\code{\link{statistics}}}
\arguments{
  \item{x}{A sampler}
  \item{names}{Which statistics to extract? A character vector of
statistics collected during the MCMC run.}
  \item{...}{Further arguments passed to method}
}
