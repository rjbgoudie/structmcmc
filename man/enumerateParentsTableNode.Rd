\name{enumerateParentsTableNode}
\alias{enumerateParentsTableNode}
\title{Table of all possible parents of a node...}
\usage{enumerateParentsTableNode(node, numberOfNodes, maxNumberParents,
    required, banned)}
\description{Table of all possible parents of a node}
\details{Creates a matrix, with each row being a parent set.
The is is created subject to the supplied indegree restriction, and the 
the supplied \code{required} and \code{banned} restrictions.}
\value{A matrix with \code{maxNumberParents} columns.
Each row is a possible parent set for node \code{node}, accounting for 
the restrictions given by \code{required} and \code{banned}.
Entries that are \code{NA} indicate no parent. e.g. there is only one 
parent, the other entries will be \code{NA}.}
\arguments{\item{node}{A node}
\item{numberOfNodes}{Number of nodes}
\item{maxNumberParents}{The maximum indegree of the node.}
\item{required}{A numeric vector. The nodes that must be included in 
every parent set for this node.}
\item{banned}{A numeric vector. The nodes that must be excluded from 
every parent set for this node.}}
