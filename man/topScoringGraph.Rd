\name{topScoringGraph}
\alias{topScoringGraph}
\title{Get modal graph given an order.}
\usage{
  topScoringGraph(order, numberOfNodes, nodesSeq,
  scoresParents, parentsTables, allRows, rowsThatContain)
}
\arguments{
  \item{order}{A vector length \code{numberOfNodes}, giving
  a permuation of \code{1:numberOfNodes}.}

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
  Returns the modal network. A \code{currentNetwork}
  object.
}
\description{
  Get modal graph given an order.
}

