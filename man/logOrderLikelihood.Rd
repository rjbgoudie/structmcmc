\name{logOrderLikelihood}
\alias{logOrderLikelihood}
\title{Log likelihood of an order...}
\usage{logOrderLikelihood(order, numberOfNodes, nodesSeq,
    scoresParents, parentsTables, allRows, rowsThatContain)}

\description{
  Log likelihood of an order
}

\details{
  Evaluates the log likelihood of an order, using the
  efficient decompostion used by Freidman and Koller (2003)
  (equation 8).
}
\value{Returns the sampled network. A \code{currentNetwork} object.}
\seealso{\code{\link{BNGibbsSampler}}, \code{\link{samplePair}}}
\references{Friedman, N., & Koller, D. (2003). \emph{Being Bayesian about Network
Structure. A Bayesian Approach to Structure Discovery in Bayesian
Networks}. Machine Learning, 50, 95-125.
\url{http://dx.doi.org/10.1023/A:1020249912095}}
\arguments{
  \item{order}{A vector length \code{numberOfNodes}, giving a permuation
of \code{1:numberOfNodes}.}
  \item{numberOfNodes}{The number of nodes in the network. A numeric vector 
of length 1.}
  \item{nodesSeq}{The vector 1:nNodes(currentNetwork). (Supplied as an 
argument for possible speed gain)}
  \item{scoresParents}{A list of the form returned by 
\code{scoreParentsTable()}}
  \item{parentsTables}{A list of tables of the form returned by 
\code{enumerateParentsTable()}}
  \item{allRows}{The vector 1:nrow(parentsTables). (Supplied as an 
argument for possible speed gain)}
  \item{rowsThatContain}{A list of the form created by 
\code{getRowsThatContain()}}
}
