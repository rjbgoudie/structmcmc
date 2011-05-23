\name{eachChangesChoicesForRequired}
\alias{eachChangesChoicesForRequired}
\title{Each changes choices for required.}
\usage{eachChangesChoicesForRequired(i, possibleParents,
    currentNetwork, change, descendantsList)}

\description{
  Each changes choices for required.
}
\value{A table(?) of options}
\arguments{
  \item{i}{The node}
  \item{possibleParents}{The output of \code{\link{getPossibleParents}}}
  \item{currentNetwork}{A list, containing in the first position the
starting \code{bn}, and in the second position the routes matrix for
that BN.}
  \item{change}{A numeric vector, containing the nodes whose parents
are to be changed.}
  \item{descendantsList}{The output of \code{\link{allDescendants}}}
}
