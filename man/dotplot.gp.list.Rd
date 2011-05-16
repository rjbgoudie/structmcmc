\name{dotplot.gp.list}
\alias{dotplot.gp.list}
\title{Dotplot of posterior graph probabilities.}
\usage{\method{dotplot}{gp.list} (gplist, subsetBy="Exact", head=30, ...)}
\description{Dotplot of posterior graph probabilities.}
\details{method description}
\seealso{\code{\link{xyplot.gp.list}}, \code{\link{dotplot.gp}}}
\arguments{\item{gplist}{...}
\item{subsetBy}{...}
\item{head}{...}
\item{...}{Further arguments passed to method}}
\examples{x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
x <- data.frame(x1 = x1, x2 = x2)

exact <- posterior(data = x, "exact")
mcmc <- posterior(data = x, "mc3", nSamples = 500, nBurnin = 100)

mygp1 <- gp(exact)
mygp2 <- gp(mcmc)
if (require(lattice)){
dotplot(gp.list(Exact = mygp1, MCMC = mygp2))
}}
