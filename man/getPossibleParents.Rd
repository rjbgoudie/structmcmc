\name{getPossibleParents}
\alias{getPossibleParents}
\title{Get possible parents.}
\usage{getPossibleParents(bn, nonDescendantsList, descendantsList,
    numberOfNodes, change)}
\description{Get possible parents.}
\details{Give a sub-bn, get the possible parents of each node.}
\value{A list of possible parents of each node}
\arguments{\item{bn}{A sub-bn}
\item{nonDescendantsList}{The output of \code{\link{nonDescendantsList}}}
\item{descendantsList}{The output of \code{\link{descendantsList}}}
\item{numberOfNodes}{The number of nodes of the bn}
\item{change}{A numeric vector, containing the nodes whose parents
are to be changed.}}
