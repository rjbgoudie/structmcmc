\name{parentalToCPDAG}
\alias{parentalToCPDAG}
\title{Convert 'parental list' to CPDAGs.}
\usage{parentalToCPDAG(x, verbose=T)}

\description{
  Convert 'parental list' to CPDAGs.
}

\details{
  A wrapper for \code{\link[parental]{as.cpdag}} that
  routes around the issue that the parental.list is not of
  class \code{bn.list} etc.
}
\seealso{\code{\link[parental]{as.cpdag}}, the function for which this
function is a wrapper.}
\value{A \code{parental.list} of CPDAGs.}
\arguments{
  \item{x}{A \code{parental.list}}
  \item{verbose}{A logical. Should a progress bar be displayed?}
}
