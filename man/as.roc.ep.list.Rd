\name{as.roc.ep.list}
\alias{as.roc.ep.list}
\title{Prepare an list of edge probability matrix for a ROC plot.}
\usage{\method{as.roc}{ep.list} (x, true, thresholds, labels, verbose, ...)}
\description{Prepare an list of edge probability matrix for a ROC plot.}
\details{Compares a list of edge probability matrices \code{x} to a true graph
\code{true}, and returns the number of true and false positives.}
\seealso{For individual edge probability matrices see
\code{\link{as.roc.ep}}. For graphs, \code{\link{as.roc.parental}} and
\code{\link{as.roc.parental.list}}.}
\value{A data frame, with columns \code{estimate}, \code{tp}, and
\code{fp}. The first contains the supplied label; the latter two
contain the number of true and false positives respectively. Rows
correspond to the supplied thresholds, for each element of the supplied
\code{ep.list}.}
\arguments{\item{x}{A list of matrices of class \code{ep} of edge probabilities}
\item{true}{A \code{bn}. The true graph.}
\item{thresholds}{A numeric vector of thresholds.}
\item{labels}{A vector of labels}
\item{verbose}{A logical. Should progress should be displayed?}
\item{...}{Further arguments (unused)}}
