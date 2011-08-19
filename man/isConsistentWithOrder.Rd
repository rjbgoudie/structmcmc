\name{isConsistentWithOrder}
\alias{isConsistentWithOrder}
\title{Is a BN consistent with an order?}
\usage{
  isConsistentWithOrder(x, order)
}
\arguments{
  \item{x}{A \code{bn}.}

  \item{order}{A vector length \code{numberOfNodes}, giving
  a permuation of \code{1:numberOfNodes}.}
}
\value{
  A logical indicating whether the BN is consistent with
  the supplied ordering.
}
\description{
  Tests whether a graph is consistent with a node ordering.
}

