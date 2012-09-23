\name{allBannedExceptPPNotBannedIfNotRequired}
\alias{allBannedExceptPPNotBannedIfNotRequired}
\title{Make a banned list.}
\usage{
  allBannedExceptPPNotBannedIfNotRequired(nodesSeq,
    possibleParents, change, toBanIfNotRequired, required)
}
\arguments{
  \item{nodesSeq}{A sequence \code{1:numberOfNodes}}

  \item{possibleParents}{The output of
  \code{\link{getPossibleParents}}}

  \item{change}{A numeric vector, containing the nodes
  whose parents are to be changed.}

  \item{toBanIfNotRequired}{to ban if not required}

  \item{required}{A list of required}
}
\value{
  A list
}
\description{
  All banned except possible parents that are not required
}

