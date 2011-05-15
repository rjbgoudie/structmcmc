\name{logNumMHNeighbours}
\alias{logNumMHNeighbours}
\title{(Log) Number of neighbouring networks.}
\usage{logNumMHNeighbours(routes, adjacency, constraintT,
    maxNumberParents=Inf)}
\description{(Log) Number of neighbouring networks.}
\details{Returns the number of acyclic graphs that can be formed by adding, 
removing or flipping a single edge of the current network}
\value{A numeric of length 1. The log number of neighbouring graphs.}
\arguments{\item{routes}{The routes matrix of the network}
\item{adjacency}{The adjacency matrix of the network}
\item{constraintT}{The transpose of a constraint matrix}
\item{maxNumberParents}{Integer of length 1. The maximum number of
parents of any node.}}
