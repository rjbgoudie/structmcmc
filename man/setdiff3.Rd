\name{setdiff3}
\alias{setdiff3}
\title{Fast, dangerous set difference.}
\usage{
  setdiff3(x, y)
}
\arguments{
  \item{x}{A vector, of the same mode as \code{y}.}

  \item{y}{A vector, of the same mode as \code{x}.}
}
\value{
  A vector of the same mode as the inputs.
}
\description{
  A fast, simple version of \code{setdiff}. This version
  does not handle factors.
}
\details{
  Note that no sanity checks are performed on the input.
}

