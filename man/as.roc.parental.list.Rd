\name{as.roc.parental.list}
\alias{as.roc.parental.list}
\title{Prepare a parental list for a ROC plot.}
\usage{
  \method{as.roc}{parental.list} (x, true, labels, verbose,
    ...)
}
\arguments{
  \item{x}{The object to be prepared for an ROC plot}

  \item{true}{A \code{bn}. The true graph.}

  \item{labels}{A vector of labels}

  \item{verbose}{A logical. Should progress be displayed?}

  \item{...}{Further arguments (unused)}
}
\value{
  A data frame, with columns \code{estimate}, \code{tp},
  and \code{fp}. The first contains the supplied label; the
  latter two contain the number of true and false positives
  respectively
}
\description{
  Compares a list of graphs \code{x} to a true graph
  \code{true}, and returns the number of true and false
  positives.
}
\seealso{
  For individual graphs, see \code{\link{as.roc.parental}}.
  For edge probability matrices see \code{\link{as.roc.ep}}
  and \code{\link{as.roc.ep.list}}.
  \code{\link{as.roc.parental}}
}

