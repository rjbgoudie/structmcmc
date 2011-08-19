\name{whichNum}
\alias{whichNum}
\title{Find rows satisfying equality.}
\usage{
  whichNum(col, x, i)
}
\arguments{
  \item{col}{The column of the matrix to use. A numeric
  vector of length 1.}

  \item{x}{A matrix.}

  \item{i}{The value to check for}
}
\value{
  A numeric vector, of the
}
\description{
  Find which rows in column \code{col} of a matrix \code{x}
  are equal to \code{i}.
}
\details{
  Note that no sanity checks are performed on the input.
}

