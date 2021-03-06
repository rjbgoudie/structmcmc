\name{as.roc}
\alias{as.roc}
\title{Prepare data for plotting in a ROC plot.}
\usage{
  as.roc(x, ...)
}
\arguments{
  \item{x}{The object to be prepared for an ROC plot}

  \item{...}{Further arguments passed to method}
}
\value{
  A data frame, with columns \code{estimate}, \code{tp},
  and \code{fp}. The first contains the supplied label; the
  latter two contain the number of true and false positives
  respectively.
}
\description{
  Converts MCMC samples, edge probability matrices etc to a
  data frame ready for plotting as ROC curves.
}
\seealso{
  For actual graphs, see \code{\link{as.roc.parental}} and
  \code{\link{as.roc.parental.list}}. For edge probability
  matrices see \code{\link{as.roc.ep}} and
  \code{\link{as.roc.ep.list}}.  For plotting
  \code{\link{rocplot}}. The methods defined are
  \code{\link{as.roc.parental}},
  \code{\link{as.roc.parental.list}},
  \code{\link{as.roc.ep}}, \code{\link{as.roc.ep.list}}
}

