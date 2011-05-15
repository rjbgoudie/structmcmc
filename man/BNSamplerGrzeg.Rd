\name{BNSamplerGrzeg}
\alias{BNSamplerGrzeg}
\title{Undocumented.}
\usage{BNSamplerGrzeg(data, initial, prior, return="network",
    logScoreFUN=defaultLogScoreFUN(),
    logScoreParameters=list(hyperparameters = "qi"), constraint,
    modejumping=F, verbose=F, keepTape=F,
    grzeg=list(grzegMoveProbability = 0.5))}
\description{Undocumented.}
\details{method description. Not tested.}
\seealso{\code{\link{BNSampler}}, \code{\link{BNSamplerPT}},
\code{\link{BNSamplerMJ}}, \code{\link{BNGibbsSampler}}.
Internally uses \code{\link{localPartitionFunction}}, and 
\code{\link{scoreParents}}.}
\arguments{\item{data}{...}
\item{initial}{...}
\item{prior}{...}
\item{return}{...}
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
\item{logScoreParameters}{...}
\item{constraint}{...}
\item{modejumping}{NOT IMPLEMENTED}
\item{verbose}{...}
\item{keepTape}{...}
\item{grzeg}{...}}
