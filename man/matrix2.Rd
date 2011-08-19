\name{matrix2}
\alias{matrix2}
\title{Fast, dangerous matrix generation.}
\usage{
  matrix2(data, nrow, ncol)
}
\arguments{
  \item{data}{A data vector.}

  \item{nrow}{The desired number of rows}

  \item{ncol}{The desired number of columns}
}
\value{
  A matrix of dimension \code{nrow} by \code{ncol},
  containing data \code{data}.
}
\description{
  A fast, simple version of \code{matrix}.
}
\details{
  The matrix must be supplied by column (i.e. \code{byrow =
  FALSE}) and the matrix can not have any names (i.e.
  \code{dimnames = NULL}).

  Note that no sanity checks are performed on the input.
}

