\name{numberDAGsGivenOrder}
\alias{numberDAGsGivenOrder}
\title{Number of BNs consistent with an order.}
\usage{
  numberDAGsGivenOrder(fam, order)
}
\arguments{
  \item{fam}{A \code{parental.list}}

  \item{order}{A vector length \code{numberOfNodes}, giving
  a permuation of \code{1:numberOfNodes}.}
}
\value{
  A numeric. The number of BNs consistent with the order
}
\description{
  Computes the number of BNs from a list of BNs that are
  consistent with an order
}

