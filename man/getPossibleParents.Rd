\name{getPossibleParents}
\alias{getPossibleParents}
\title{Get possible parents.}
\usage{getPossibleParents(bn, nonDescendantsList, descendantsList,
    numberOfNodes, change, maxIndegree)}

\description{
  Get possible parents.
}

\details{
  Given a sub-bn, get the possible parents of each node.

  The possible parents for a node j that changes (a "change
  node") are those that are not descendants of any other
  change node, except for any change node that is an
  immediate parents of j in bn (the "sub-bn").
}
\value{A list of possible parents of each node}
\arguments{
  \item{bn}{A sub-bn}
  \item{nonDescendantsList}{The output of \code{\link{allNonDescendants}}}
  \item{descendantsList}{The output of \code{\link{allDescendants}}}
  \item{numberOfNodes}{The number of nodes of the bn}
  \item{change}{A numeric vector, containing the nodes whose parents
are to be changed.}
  \item{maxIndegree}{Maximum indegree}
}
