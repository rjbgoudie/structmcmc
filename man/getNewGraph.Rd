\name{getNewGraph}
\alias{getNewGraph}
\title{List graphs in change node neighbourhood.}
\usage{getNewGraph(currentNetwork, change)}

\description{
  List graphs in change node neighbourhood.
}

\details{
  Given a list of nodes \code{change}, get a list of
  \code{bn}.
}
\value{A list of \code{bn}}
\arguments{
  \item{currentNetwork}{A list, containing in the first position the
starting \code{bn}, and in the second position the routes matrix for
that BN.}
  \item{change}{A numeric vector, containing the nodes whose parents
are to be changed.}
}
