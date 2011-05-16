\name{BNSamplerMJ}
\alias{BNSamplerMJ}
\title{Mode-jumping MCMC sampler for Bayesian Networks.}
\usage{BNSamplerMJ(data, initial, prior, return="network",
    logScoreFUN=logScoreMultDirFUN(),
    logScoreParameters=list(hyperparameters = "qi"), constraint,
    modejumping=F, verbose=F, keepTape=F)}
\description{Mode-jumping MCMC sampler for Bayesian Networks.}
\details{The sampler samples Bayesian Networks (ie models).}
\value{A function, which when called draws the next sample of the MCMC.}
\seealso{\code{\link{BNSampler}}, \code{\link{BNGibbsSampler}},
\code{\link{BNSamplerPT}}, \code{\link{BNSamplerGrzeg}},
\code{\link{BNSamplerBigFlips}}. Internally uses 
\code{\link{whichGraphsNotNeighbours}}}
\arguments{\item{data}{The data.}
\item{initial}{An object of class 'bn'. The starting value of the MCMC.}
\item{prior}{A function that returns the prior score of the supplied bn.}
\item{return}{Either "network" or "contingency".}
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
}}
\item{logScoreParameters}{A list of parameters that are passed to
logScoreFUN.}
\item{constraint}{A matrix of dimension ncol(data) x ncol(data) giving
constraints to the sample space. The (i, j) element is:
\describe{
\item{1}{if the edge i -> j is required}
\item{-1}{if the edge i -> is excluded.}
\item{0}{if the edge i -> j is not constrained.}
}
The diagonal of constraint must be all 0.}
\item{modejumping}{Either a logical of length 1, or a list. When no
mode-jumping is desired, use modejumping = F. For mode-jumping, use a 
list with the following components:
\describe{
\item{modes}{A bn.list of modes.}
\item{modesLogScores}{Optionally, a numeric vector, containing the 
logScores of the modes.}
\item{checkModesAcyclic}{A logical of length 1}
\item{modesPreFiltered}{A logical of length 1}
\item{modeJumpingProbability}{A numeric of length 1. Default 0.25}
\item{dontCheckModesValid}{A logical of length 1}
}}
\item{verbose}{A logical of length 1, indicating whether verbose
output should be printed.}
\item{keepTape}{A logical of length 1, indicating whether a full log
(\code{tape}) of the MCMC sampler should be kept. Enabling this option 
can be very memory-intensive.}}
