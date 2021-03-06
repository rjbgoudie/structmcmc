\name{validStartEnd}
\alias{validStartEnd}
\title{Check validity of start and end.}
\usage{
  validStartEnd(start, end, length)
}
\arguments{
  \item{start}{The start value to check}

  \item{end}{The end value to check}

  \item{length}{A numeric of length 1. The length of the
  MCMC run that start and end should window.}
}
\value{
  A logical of length 1 indicating if start and end are
  valid
}
\description{
  Checks if start and end are valid as start and end points
  for taking a window (a subset, e.g. to remove burn-in) in
  a MCMC run
}

