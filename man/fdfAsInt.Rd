\name{fdfAsInt}
\alias{fdfAsInt}
\title{Convert a data frame of factors to integers.}
\usage{fdfAsInt(x, useLevelNames=T)}

\description{
  Convert a data frame of factors to integers.
}

\details{
  Converts a \code{data.frame} that consists columns of
  factors into a \code{data.frame} consisting of integers
}
\value{The data.frame with columns converted to integers}
\seealso{\code{\link{intAsFDF}}}
\arguments{
  \item{x}{an object of class \code{data.frame}}
  \item{useLevelNames}{logical, indicating whether the labels of the 
levels should be converted to integers. This only makes sense if the 
levels are integers stored as characters. e.g. factor(c("3", "2", "3"))}
}
