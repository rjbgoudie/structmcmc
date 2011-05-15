\name{as.roc.parental}
\alias{as.roc.parental}
\title{Prepare a parental for a ROC plot.}
\usage{\method{as.roc}{parental} (x, true, label, ...)}
\description{Prepare a parental for a ROC plot.}
\details{Compares two graphs, one of which is the true graph, and computes the
number of true and false positives.}
\seealso{For lists of graphs, see \code{\link{as.roc.parental.list}}. For
edge probability matrices see \code{\link{as.roc.ep}} and
\code{\link{as.roc.ep.list}}.}
\value{A 1-row data frame, with columns \code{estimate}, \code{tp}, and
\code{fp}. The first contains the supplied label; the latter two
contain the number of true and false positives respectively}
\seealso{\code{\link{as.roc.parental.list}}}
\arguments{\item{x}{A \code{parental} graph to compare to \code{true}}
\item{true}{A \code{bn}. The true graph.}
\item{label}{A label}
\item{...}{Further arguments (unused)}}
