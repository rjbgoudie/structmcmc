\name{scoreParentsTableNode}
\alias{scoreParentsTableNode}
\title{Score a node-level parentsTable...}
\usage{scoreParentsTableNode(node, parentsTable, logScoreOfflineFUN,
    logScoreParameters, prior)}
\description{Score a node-level parentsTable}
\details{Computes the local scores of all the Bayesian Networks, with node 
\code{node} set to parent sets corresponding to each row of a single 
component of a \code{parentsTable}.}
\value{A numeric vector of scores.}
\arguments{\item{node}{A node. A numeric vector of length 1.}
\item{parentsTable}{A component of a \code{parentsTable}, of the form 
created by \code{enumerateParentsTable()}.}
\item{logScoreOfflineFUN}{A function that computes the logScore
of a Bayesian Network.}
\item{logScoreParameters}{A list of parameters that are passed to
\code{logScoreFUN}.}
\item{prior}{function that returns the prior score of the
supplied bn.}}
