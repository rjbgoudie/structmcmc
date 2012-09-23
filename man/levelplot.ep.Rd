\name{levelplot.ep}
\alias{levelplot.ep}
\title{Levelplot of posterior edge probabilities.}
\usage{
  \method{levelplot}{ep} (ep)
}
\arguments{
  \item{ep}{An object of class \code{ep}. A matrix of edge
  probabilities.}
}
\description{
  Plot a \code{\link[lattice]{levelplot}} displaying the
  edge probabilities.
}
\examples{
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
x <- data.frame(x1 = x1, x2 = x2)

exact <- posterior(data = x, "exact")

myep <- ep(exact)
if (require(lattice)){
  levelplot(myep)
}
}
\seealso{
  \code{\link{levelplot.ep.list}}, \code{\link{dotplot.ep}}
}

