\name{scoreParentsTable}
\alias{scoreParentsTable}
\title{Score a node-level parentsTable...}
\usage{scoreParentsTable(parentsTables, logScoreOfflineFUN,
    logScoreParameters, prior)}
\description{Score a node-level parentsTable}
\details{Computes the scores of all the Bayesian Networks, with parent sets 
corresponding to each row of a single component of a \code{parentsTable}.}
\value{List of numeric vectors of scores.}
\arguments{\item{parentsTables}{A component of a \code{parentsTables}, of the form 
created by \code{enumerateParentsTable()}.}
\item{logScoreOfflineFUN}{A function that computes the logScore
of a Bayesian Network.}
\item{logScoreParameters}{A list of parameters that are passed to
\code{logScoreFUN}.}
\item{prior}{function that returns the prior score of the
supplied bn.}}
