\name{requireSomethingFromEachParent}
\alias{requireSomethingFromEachParent}
\title{Require something from each parent.}
\usage{
  requireSomethingFromEachParent(numberOfNodes, nodesSeq,
    possibleParents, change, descendantsList,
    currentNetwork, maxIndegree)
}
\arguments{
  \item{numberOfNodes}{The number of nodes}

  \item{nodesSeq}{A sequence \code{1:numberOfNodes}}

  \item{possibleParents}{The output of
  \code{\link{getPossibleParents}}}

  \item{change}{A numeric vector, containing the nodes
  whose parents are to be changed.}

  \item{descendantsList}{The output of
  \code{\link{allDescendants}}}

  \item{currentNetwork}{A list, containing in the first
  position the starting \code{bn}, and in the second
  position the routes matrix for that BN.}

  \item{maxIndegree}{Maximum indegree}
}
\value{
  A list
}
\description{
  Require something from each parent.
}

