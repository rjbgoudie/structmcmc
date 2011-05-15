\name{xyplot.gp}
\alias{xyplot.gp}
\title{Scatterplot of posterior graph probabilities.}
\usage{\method{xyplot}{gp} (gp, head=30, ...)}
\description{Scatterplot of posterior graph probabilities.}
\details{method description}
\seealso{\code{\link{dotplot.gp}}}
\arguments{\item{gp}{...}
\item{head}{...}
\item{...}{Further arguments passed to method}}
\examples{x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
x <- data.frame(x1 = x1, x2 = x2)

exact <- posterior(data = x, "exact")

mygp <- gp(exact)
if (require(lattice)){
xyplot(mygp)
}}
