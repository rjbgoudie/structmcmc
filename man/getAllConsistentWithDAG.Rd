\name{getAllConsistentWithDAG}
\alias{getAllConsistentWithDAG}
\title{Get all consistent with DAG}
\usage{
  getAllConsistentWithDAG(bn, currentNetwork,
  numberOfNodes, nodesSeq, change, maxIndegree)
}
\arguments{
  \item{bn}{A sub-bn}

  \item{currentNetwork}{A list, containing in the first
  position the starting \code{bn}, and in the second
  position the routes matrix for that BN.}

  \item{numberOfNodes}{The number of nodes of the bn}

  \item{nodesSeq}{A sequence \code{1:numberOfNodes}}

  \item{change}{A numeric vector, containing the nodes
  whose parents are to be changed.}

  \item{maxIndegree}{Maximum indegree}
}
\value{
  A list of possible parents of each node
}
\description{
  Get all consistent with DAG
}

