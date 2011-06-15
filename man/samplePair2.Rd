\name{samplePair2}
\alias{samplePair2}
\title{Sample the parents of a pair of nodes (Gibbs sampler) v2.}
\usage{samplePair2(currentNetwork, numberOfNodes, nodesSeq,
    scoresParents, parentsTables, allRows, rowsThatContain,
    logScoreFUN, logScoreParameters)}

\description{
  Sample the parents of a pair of nodes (Gibbs sampler) v2.
}

\details{
  Sample from posterior distribution on graph, conditional
  on all the edges, except for those corresponding to the
  parents sets of two nodes.
}
\value{Returns the sampled network. A \code{currentNetwork} object.}
\seealso{\code{\link{BNGibbsSampler}}, \code{\link{sampleNode}}}
\arguments{
  \item{currentNetwork}{A \code{currentNetwork} object}
  \item{numberOfNodes}{The number of nodes in the network. A numeric vector
of length 1.}
  \item{nodesSeq}{The vector 1:nNodes(currentNetwork). (Supplied as an
argument for possible speed gain)}
  \item{scoresParents}{A list of the form returned by
\code{scoreParentsTable()}}
  \item{parentsTables}{A list of tables of the form returned by
\code{enumerateParentsTable()}}
  \item{allRows}{The vector 1:nrow(parentsTables). (Supplied as an
argument for possible speed gain)}
  \item{rowsThatContain}{A list of the form created by
\code{getRowsThatContain()}}
}
