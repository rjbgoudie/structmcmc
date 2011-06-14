\name{orderFlipOperator}
\alias{orderFlipOperator}
\title{Draw order proposal using flip-operator...}
\usage{orderFlipOperator(order, nodesSeq)}

\description{
  Draw order proposal using flip-operator
}

\details{
  Draws a proposal order, given a current order, using the
  single flip-operator. See Grzegorcyck & Husmeier (2008)
  eq 12.
}
\value{A new order}
\arguments{
  \item{order}{A permutation of \code{nodesSeq}}
  \item{nodesSeq}{The vector 1:nNodes(currentNetwork). (Supplied as an 
argument for possible speed gain)}
}
