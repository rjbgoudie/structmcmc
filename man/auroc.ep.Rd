\name{auroc.ep}
\alias{auroc.ep}
\title{Compute the area under an ROC curve.}
\usage{
  \method{auroc}{ep} (x, true, thresholds, label, verbose,
  ...)
}
\arguments{
  \item{x}{A matrix of class \code{ep} of edge
  probabilities}

  \item{true}{A \code{bn}. The true graph.}

  \item{thresholds}{A numeric vector of thresholds.}

  \item{label}{A label}

  \item{verbose}{A logical. Should progress should be
  displayed?}

  \item{...}{Further arguments (unused)}
}
\value{
  The area under the ROC curve.
}
\description{
  Compute the area under an ROC curve for an edge
  probability matrix.
}
\seealso{
  \code{\link{auroc}}, \code{\link{auroc.ep.list}}
}

