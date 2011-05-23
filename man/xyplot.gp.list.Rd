\name{xyplot.gp.list}
\alias{xyplot.gp.list}
\title{Scatterplot of posterior graph probablities.}
\usage{\method{xyplot}{gp.list} (gplist, head=30, scales=list(x =
    list(rot = 90)), highlight, ...)}

\description{
  Scatterplot of posterior graph probablities.
}

\details{
  method description
}
\seealso{\code{\link{dotplot.gp.list}}, \code{\link{xyplot.gp}}}
\arguments{
  \item{gplist}{...}
  \item{head}{...}
  \item{scales}{...}
  \item{highlight}{...}
  \item{...}{Further arguments passed to method}
}
\examples{x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
x <- data.frame(x1 = x1, x2 = x2)

exact <- posterior(data = x, "exact")
mcmc <- posterior(data = x, "mc3", nSamples = 500, nBurnin = 100)

mygp1 <- gp(exact)
mygp2 <- gp(mcmc)
if (require(lattice)){
xyplot(gp.list(Exact = mygp1, MCMC = mygp2))
}}
