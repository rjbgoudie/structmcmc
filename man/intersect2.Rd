\name{intersect2}
\alias{intersect2}
\title{Fast, dangerous set intersect.}
\usage{
  intersect2(x, y, nmax = NA)
}
\arguments{
  \item{x}{A vector, of the same mode as \code{y}.}

  \item{y}{A vector, of the same mode as \code{x}.}
}
\value{
  A vector of the same mode as the inputs.
}
\description{
  A fast, simple version of \code{intersect}. This version
  does not handle factors.
}
\details{
  Note that no sanity checks are performed on the input.
}

