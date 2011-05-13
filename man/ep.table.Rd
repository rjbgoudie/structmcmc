\name{ep.table}
\alias{ep.table}
\title{Computes the edge probabilities implied by a table.}
\usage{\method{ep}{table} (x, verbose=F, ...)}
\description{Computes the edge probabilities implied by a table.}
\details{Computes the edge probabilities implied by a table.}
\value{A matrix of class 'ep' with entry (i,j) containing the probability
of an edge from node i --> j}
\arguments{\item{x}{An object of class 'table'. The table should have names that
can be converted into objects of class parental using
\code{as.parental.character}. ie they should be serializations of
parental objects.}
\item{verbose}{A logical indicating whether progress should be shown.}
\item{...}{Further arguments passed to
\code{\link{ep.parental.contingency}}}}
