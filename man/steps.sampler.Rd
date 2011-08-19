\name{steps.sampler}
\alias{steps.sampler}
\title{Number of samples drawn.}
\usage{
  \method{length}{sampler} (x, ...)
}
\arguments{
  \item{x}{A sampler}

  \item{...}{Further arguments, currently unused}
}
\value{
  The number of samples (steps) that have been drawn.
}
\description{
  Returns the number of samples (MCMC steps) drawn in the
  supplied sampler.
}

