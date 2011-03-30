\name{ep.parental.contingency}
\alias{ep.parental.contingency}
\title{Computes the edge probabilities...}
\usage{\method{ep}{parental.contingency} (x, FUN, verbose=F, ...)}
\description{Computes the edge probabilities}
\details{Computes the edge probabilities from a \code{parental.contingency},
which is a list, whose first component is \code{parental.list}, and
whose second component is a contingency table of counts of the
corresponding element of the \code{parental.list}. The contingency
table must simply be a numeric vector. The first component must be named
\code{parental.list}, and the second \code{contingency}.}
\value{A matrix of class 'ep' with entry (i,j) containing the
probability of an edge from node \code{i} to \code{j}}
\arguments{\item{x}{An object of class \code{parental.contingency}, of the form
desribed in the Details section.}
\item{FUN}{A function to be applied to the samples before the edge
probabilities are computed. The function should accept a vector
of objects of class 'parental.list' and return the transformed
parental objects as a single 'parental.list'. If the function
does not return a 'parental.list', an error will be thrown. Can be
supplied as a function, symbol or character string (see
\code{\link[base]{match.fun}})}
\item{verbose}{A logical indicating whether progress should be shown.}
\item{...}{Further arguments (unused)}}
