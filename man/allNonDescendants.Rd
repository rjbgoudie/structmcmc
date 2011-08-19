\name{allNonDescendants}
\alias{allNonDescendants}
\title{Get non-descendants of all nodes.}
\usage{
  allNonDescendants(currentNetwork)
}
\arguments{
  \item{currentNetwork}{A list, containing in the first
  position the starting \code{bn}, and in the second
  position the routes matrix for that BN.}
}
\value{
  A list, each componentof which contains a list of the
  non-descendants of the corresponding node of the supplied
  BN.
}
\description{
  Get non-descendants of all nodes.
}

