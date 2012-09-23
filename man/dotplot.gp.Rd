\name{dotplot.gp}
\alias{dotplot.gp}
\title{Dotplot of posterior graph probabilities.}
\usage{
  \method{dotplot}{gp} (gp, head = 30, ...)
}
\arguments{
  \item{gp}{...}

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

mygp <- gp(exact)
if (require(lattice)){
  dotplot(mygp)
}
}
\seealso{
  \code{\link{xyplot.gp}}
}

