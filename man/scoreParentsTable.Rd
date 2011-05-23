\name{scoreParentsTable}
\alias{scoreParentsTable}
\title{Score a node-level parentsTable.}
\usage{scoreParentsTable(parentsTables, logScoreLocalFUN,
    logScoreParameters, prior, verbose=F)}

\description{
  Score a node-level parentsTable.
}

\details{
  Computes the scores of all the Bayesian Networks, with
  parent sets corresponding to each row of a single
  component of a \code{parentsTable}.
}
\value{List of numeric vectors of scores.}
\seealso{\code{\link{scoreParentsTableNode}},
\code{\link{enumerateParentsTable}}, \code{\link{whichParentSetRows}},
\code{\link{getRowsThatContain}}}
\arguments{
  \item{parentsTables}{A component of a \code{parentsTables}, of the form 
created by \code{enumerateParentsTable()}.}
  \item{logScoreLocalFUN}{A function that computes the local logScore 
of a Bayesian Network.}
  \item{logScoreParameters}{A list of parameters that are passed to
\code{logScoreFUN}.}
  \item{prior}{function that returns the prior score of the
supplied bn.}
  \item{verbose}{A logical of length 1, indicating whether verbose output
should be printed.}
}
