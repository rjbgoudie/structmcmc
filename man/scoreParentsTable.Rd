\name{scoreParentsTable}
\alias{scoreParentsTable}
\title{Score a node-level parentsTable.}
\usage{
  scoreParentsTable(parentsTables, logScoreLocalFUN,
    logScoreParameters, prior, verbose = F)
}
\arguments{
  \item{parentsTables}{A component of a
  \code{parentsTables}, of the form created by
  \code{enumerateParentsTable()}.}

  \item{logScoreLocalFUN}{A function that computes the
  local logScore of a Bayesian Network.}

  \item{logScoreParameters}{A list of parameters that are
  passed to \code{logScoreFUN}.}

  \item{loalPriors}{A list of functions of the same length
  as the number of nodes that returns the local prior score
  of the corresponding node, given a numeric vector of
  parents.}

  \item{verbose}{A logical of length 1, indicating whether
  verbose output should be printed.}
}
\value{
  List of numeric vectors of scores.
}
\description{
  Computes the scores of all the Bayesian Networks, with
  parent sets corresponding to each row of a single
  component of a \code{parentsTable}.
}
\seealso{
  \code{\link{scoreParentsTableNode}},
  \code{\link{enumerateParentsTable}},
  \code{\link{whichParentSetRows}},
  \code{\link{getRowsThatContain}}
}

