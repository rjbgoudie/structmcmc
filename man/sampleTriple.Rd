\name{sampleTriple}
\alias{sampleTriple}
\title{Sample the parents of a triple of nodes (Gibbs sampler).}
\usage{sampleTriple(currentNetwork, numberOfNodes, nodesSeq, scoresParents,
    parentsTables, allRows, rowsThatContain, logScoreFUN,
    logScoreParameters)}
\description{Sample the parents of a triple of nodes (Gibbs sampler).}
\details{Sample from posterior distribution on graph, conditional on
all the edges, except for those corresponding to the parents sets of 
two nodes.}
\value{Returns the sampled network. A \code{currentNetwork} object.}
\seealso{\code{\link{BNGibbsSampler}}, \code{\link{sampleNode}}}
\arguments{\item{currentNetwork}{A \code{currentNetwork} object}
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
\item{logScoreFUN}{A list of four elements:
\describe{
\item{offline}{A function that computes the logScore of a Bayesian 
Network}
\item{online}{A function that incrementally computes the logScore of a 
Bayesian Network}
\item{local}{A function that computes the local logScore of a 
Bayesian Network}
\item{prepare}{A function that prepares the data, and any further 
pre-computation required by the logScore functions.}
}
For Multinomial-Dirichlet models, \code{\link{logScoreMultDirFUN}} 
returns the appropriate list; for Normal models with Zellner g-priors,
\code{\link{logScoreZellnerFUN}} returns the appropriate list.}
\item{logScoreParameters}{A list of parameters that are passed to
logScoreFUN.}}
