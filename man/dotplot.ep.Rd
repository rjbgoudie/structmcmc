\name{dotplot.ep}
\alias{dotplot.ep}
\title{Dotplot of posterior edge probabilities.}
\usage{
  \method{dotplot}{ep} (ep, head = 30, ...)
}
\arguments{
  \item{ep}{...}

  \item{head}{...}

  \item{...}{Further arguments passed to method}
}
\description{
  method description
}
\examples{
x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
x <- data.frame(x1 = x1, x2 = x2)

exact <- posterior(data = x, "exact")

myep <- ep(exact)
if (require(lattice)){
dotplot(myep)
}
}
\seealso{
  \code{\link{levelplot.ep}}, \code{\link{dotplot.ep.list}}
}

