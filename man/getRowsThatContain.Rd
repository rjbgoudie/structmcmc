\name{getRowsThatContain}
\alias{getRowsThatContain}
\title{Create lookup table for parentsTable.}
\usage{
  getRowsThatContain(numberOfNodes, parentsTables,
  maxNumberParents)
}
\arguments{
  \item{numberOfNodes}{The number of nodes. A numeric
  vector.}

  \item{parentsTables}{A list of tables of the form
  returned by \code{enumerateParentsTable()}}

  \item{maxNumberParents}{The maximum indegree for each
  node. A numeric vector of length 1.}
}
\value{
  A list of length \code{numberOfNodes}. Each component of
  this list is a list of length \code{numberOfNodes}, the
  i^{th} component of which is a numeric vector containing
  the indicies of \code{parentsTables[[node]]} that
  correspond to parent sets including \code{i}.
}
\description{
  Creates a list that allows quick lookup of a
  parentsTable. This is needed for
  \code{whichParentSetRows}.
}
\seealso{
  \code{\link{enumerateParentsTable}},
  \code{\link{scoreParentsTable}}
}

