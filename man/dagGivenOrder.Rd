\name{dagGivenOrder}
\alias{dagGivenOrder}
\title{Sample a DAG given an order (weighted)}
\usage{
  dagGivenOrder(order, numberOfNodes, nodesSeq,
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
  Returns the sampled network. A \code{currentNetwork}
  object.
}
\description{
  Sample a DAG given an order, using the parent weights.
  See Grzegorcyck & Husmeier (2008), eq 15.
}

