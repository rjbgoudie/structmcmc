\name{samplePair}
\alias{samplePair}
\title{Sample the parents of a pair of nodes...}
\usage{samplePair(currentNetwork, numberOfNodes, nodesSeq, scoresParents,
    parentsTables, allRows, rowsThatContain)}
\description{Sample the parents of a pair of nodes}
\details{Sample from posterior distribution on graph, conditional on
all the edges, except for those corresponding to the parents sets of 
two nodes.}
\value{Returns the sampled network. A \code{currentNetwork} object.}
\arguments{\item{currentNetwork}{A \code{currentNetwork} object}
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
\code{getRowsThatContain()}}}
