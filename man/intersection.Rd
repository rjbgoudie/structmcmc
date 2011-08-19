\name{intersection}
\alias{intersection}
\title{Multiple intersection function}
\usage{
  intersection(x, y, ...)
}
\arguments{
  \item{x}{A vector containing a sequence of items
  (conceptually with no duplicated values).}

  \item{y}{A vector of the same \code{\link{mode}} as
  \code{x}.}

  \item{...}{Further vectors of the same mode as \code{x}.}
}
\value{
  A vector of the same \code{\link{mode}} as \code{x} and
  \code{y}.
}
\description{
  Performs intersection across all the supplied arguments.
}
\examples{
x <- c(1, 2, 4, 5)
y <- c(2, 3, 4)
intersection(x, y)
z <- c(5, 4, 3)
intersection(x, y, z)
}

