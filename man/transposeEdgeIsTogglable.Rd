\name{transposeEdgeIsTogglable}
\alias{transposeEdgeIsTogglable}
\title{Find togglable edges.}
\usage{
  transposeEdgeIsTogglable(routes, adjacency, constraintT,
    maxNumberParents = Inf)
}
\arguments{
  \item{routes}{The routes matrix of the network}

  \item{adjacency}{The adjacency matrix of the network}

  \item{constraintT}{The transpose of a constraint matrix}

  \item{maxNumberParents}{Integer of length 1. The maximum
  number of parents of any node.}
}
\value{
  A logical matrix of the same dimension as the supplied
  matrices, with entries indicating whether the
  corresponding edge can be added or removed without
  introducing a cycle. NOTE THIS IS TRANSPOSE OF EXPECTED
}
\description{
  Finds "edge-locations" in the graph that can be added or
  removed without introducing a cycle into the graph
}
\seealso{
  \code{\link{BNSampler}},
  \code{\link{transposeEdgeIsAddable}},
  \code{\link{transposeEdgeIsRemovable}},
  \code{\link{edgeIsFlippable}}
}

