\name{queryFinding}
\alias{queryFinding}
\title{Query an independence network, given a finding.}
\usage{
  queryFinding(x, node, nodeWithState, nodeState, FUN, ...)
}
\arguments{
  \item{x}{A grain}

  \item{node}{Character vector of length 1. The node whose
  marginal distribution is required.}

  \item{nodeWithState}{Character vector of length 1. The
  node whose state is to be conditioned upon.}

  \item{nodeState}{The state to set \code{nodeWithState}
  to. A character vetor of length 1.}

  \item{FUN}{a function that is applied to the output of
  \code{\link[gRain]{querygrain}}. The first argument
  should accept the result of
  \code{\link[gRain]{querygrain}} FOR A PARTICULAR NODE.
  (Usaully \code{querygrain} returns a list, since it may
  provide the distribution of a number of nodes at once.
  Here, we remove that list.)  Any other arguments are
  passed \code{...}}

  \item{...}{Passed to \code{FUN}}
}
\value{
  The result of \code{FUN}
}
\description{
  Query a finding, and apply \code{FUN}. Just a wrapper
  around \code{\link[gRain]{setFinding}} and
  \code{\link[gRain]{querygrain}}.
}
\seealso{
  \code{\link{as.grain}},
  \code{\link{marginalGivenOthers}},
  \code{\link{marginalGivenIntervention}}
}

