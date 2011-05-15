\name{whichGraphsNotNeighbours}
\alias{whichGraphsNotNeighbours}
\title{Identify neighbouring graphs.}
\usage{whichGraphsNotNeighbours(x, logScores, verbose=F, head)}
\description{Identify neighbouring graphs.}
\details{Finds which of the supplied graphs are not neighbours. Neighbouring
graphs are removed according to a heuristic involving scores.}
\value{A numeric vector, giving the indicies \code{x} which are 
not neighbours.}
\arguments{\item{x}{A list of graphs}
\item{logScores}{A numeric vector containing the log score of the 
corresponding graph in \code{x}}
\item{verbose}{Should verbose output be shown. A logical of length 1.}
\item{head}{The rank of the worst graph to return.}}
