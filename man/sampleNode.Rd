\name{sampleNode}
\alias{sampleNode}
\title{Sample the parents of a single node (Gibbs sampler).}
\usage{
  sampleNode(currentNetwork, numberOfNodes, nodesSeq,
  scoresParents, parentsTables, allRows, rowsThatContain)
}
\arguments{
  \item{currentNetwork}{A \code{currentNetwork} object}

  \item{numberOfNodes}{The number of nodes in the network.
  A numeric vector of length 1.}

  \item{nodesSeq}{The vector 1:nNodes(currentNetwork).
  (Supplied as an argument for possible speed gain)}

  \item{scoresParents}{A list of the form returned by
  \code{scoreParentsTable()}}

  \item{parentsTables}{A list of tables of the form
  returned by \code{enumerateParentsTable()}}

  \item{allRows}{The vector 1:nrow(parentsTables).
  (Supplied as an argument for possible speed gain)}

  \item{rowsThatContain}{A list of the form created by
  \code{getRowsThatContain()}}
}
\value{
  Returns the sampled network. A \code{currentNetwork}
  object.
}
\description{
  Sample from posterior distribution on graph, conditional
  on all the edges, except those that go into node
  \code{node}.
}
\seealso{
  \code{\link{BNGibbsSampler}}, \code{\link{samplePair}}
}

