\name{edgeIsFlippable}
\alias{edgeIsFlippable}
\title{Find flippable edges.}
\usage{
  edgeIsFlippable(routes, adjacency, constraintT,
    maxNumberParents)
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
  corresponding edge can be flipped without introducing a
  cycle. NOTE THIS IS TRANSPOSE OF EXPECTED
}
\description{
  Finds edges in the graph whose direction can be reversed
  ("flipped") without introducing a cycle into the graph.
}
\seealso{
  \code{\link{BNSampler}},
  \code{\link{transposeEdgeIsAddable}},
  \code{\link{transposeEdgeIsRemovable}},
  \code{\link{transposeEdgeIsTogglable}}
}

