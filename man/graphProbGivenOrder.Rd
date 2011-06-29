\name{graphProbGivenOrder}
\alias{graphProbGivenOrder}
\title{Probability of a graph given an order.}
\usage{graphProbGivenOrder(x, order, numberOfNodes, nodesSeq,
    scoresParents, parentsTables, allRows, rowsThatContain)}

\description{
  Probability of a graph given an order.
}
\value{The probabilty of that graph given the order.}
\arguments{
  \item{x}{A BN}
  \item{order}{A vector length \code{numberOfNodes}, giving a permuation
of \code{1:numberOfNodes}.}
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
