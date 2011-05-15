\name{rocplot}
\alias{rocplot}
\title{Plot an ROC curve.}
\usage{rocplot(true, maps, bnpmls, eps, thresholds=seq(from = 0, to = 1, by =
    0.01), use.cpdags=F, verbose=F, ...)}
\description{Plot an ROC curve.}
\details{Print a ROC curve given a variety of estimation methods.


An alternative is given by:
Werhli et al. Comparative evaluation of reverse engineering gene
regulatory networks with relevance networks, graphical gaussian models
and bayesian networks. Bioinformatics (2006) vol. 22 (20) pp. 2523}
\value{A lattice object, containing the ROC plot}
\seealso{\code{\link{as.roc}}}
\arguments{\item{true}{An object of class bn giving the true graph. This is
converted to a CPDAG before comparision.}
\item{maps}{An object of class \code{bn.list}. An estimate of the true
graph, given by the MAP estimate given by the MCMC. This is internally
converted to CPDAG.}
\item{bnpmls}{An object of class \code{bnpostmcmc.list}. An object
encapsulating a number of MCMC runs.}
\item{eps}{An object of class \code{ep.list}.}
\item{thresholds}{An object of class numeric, giving the thresholds at
which to plot the graphs given by MCMC and PC-MCMC estimates}
\item{use.cpdags}{A logical of length 1 indicating whether DAG should be 
converted to CPDAGs before computing the ROC curves.}
\item{verbose}{A logical indicating whether to show progress bars etc.}
\item{...}{Further arguments}}
\examples{x1 <- factor(c(1, 1, 0, 1, 1, 1))
x2 <- factor(c(0, 1, 0, 1, 0, 1))
x3 <- factor(c(0, 1, 0, 0, 0, 0))
x <- data.frame(x1 = x1, x2 = x2, x3 = x3)

exact <- posterior(x, "exact")
eppost <- ep(exact)
rocplot(true = bn(integer(0), c(1,3), integer(0)),
eps = ep.list(eppost))

mcmc <- posterior(x, "mc3", nSamples = 1000, nBurnin = 100)
rocplot(true = bn(integer(0), c(1,3), integer(0)),
eps = ep.list(eppost),
bnpmls = bnpostmcmc.list(mcmc))}
