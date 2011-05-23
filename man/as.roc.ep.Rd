\name{as.roc.ep}
\alias{as.roc.ep}
\title{Prepare an edge probability matrix for a ROC plot.}
\usage{\method{as.roc}{ep} (x, true, thresholds, label, verbose,
    ...)}

\description{
  Prepare an edge probability matrix for a ROC plot.
}

\details{
  Compares an edge probability matrices \code{x} to a true
  graph \code{true}, and returns the number of true and
  false positives.
}
\seealso{For lists of edge probability matrices, see
\code{\link{as.roc.ep.list}}. For graphs, \code{\link{as.roc.parental}}
and \code{\link{as.roc.parental.list}}. \code{\link{as.roc.ep.list}}}
\value{A data frame, with columns \code{estimate}, \code{tp}, and
\code{fp}. The first contains the supplied label; the latter two
contain the number of true and false positives respectively. Rows
correspond to the supplied thresholds.}
\arguments{
  \item{x}{A matrix of class \code{ep} of edge probabilities}
  \item{true}{A \code{bn}. The true graph.}
  \item{thresholds}{A numeric vector of thresholds.}
  \item{label}{A label}
  \item{verbose}{A logical. Should progress should be displayed?}
  \item{...}{Further arguments (unused)}
}
