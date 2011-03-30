\name{setdiff3}
\alias{setdiff3}
\title{Fast, dangerous set difference...}
\usage{setdiff3(x, y)}
\description{Fast, dangerous set difference}
\details{A fast, simple version of \code{setdiff}.
This version does not handle factors.

Note that no sanity checks are performed on the input.}
\value{A vector of the same mode as the inputs.}
\arguments{\item{x}{A vector, of the same mode as \code{y}.}
\item{y}{A vector, of the same mode as \code{x}.}}
