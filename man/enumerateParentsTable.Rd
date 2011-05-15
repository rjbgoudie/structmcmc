\name{enumerateParentsTable}
\alias{enumerateParentsTable}
\title{Table of all possible parents of a node.}
\usage{enumerateParentsTable(numberOfNodes, maxNumberParents, required,
    banned, verbose=F)}
\description{Table of all possible parents of a node.}
\details{Creates a list of tables with rows corresponding to the possible parent 
sets of each node.}
\value{A list of matrices of the form returned by 
\code{enumerateParentsTableNode}.}
\arguments{\item{numberOfNodes}{Number of nodes}
\item{maxNumberParents}{The maximum indegree of the node.}
\item{required}{A list of numeric vectors. Each component is a numeric 
vector containing the nodes that must be included in every parent set 
for this node.}
\item{banned}{A list of numeric vectors. Each component is a numeric 
vector containing the nodes that must be excluded from every parent 
set for this node.}
\item{verbose}{A logical of length 1, indicating whether verbose output
should be printed.}}
