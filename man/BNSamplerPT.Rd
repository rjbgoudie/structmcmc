\name{BNSamplerPT}
\alias{BNSamplerPT}
\title{Undocumented.}
\usage{BNSamplerPT(data, initial, prior, return="network",
    logScoreFUN=defaultLogScoreFUN(),
    logScoreParameters=list(hyperparameters = "qi"), constraint,
    modejumping=F, verbose=F, keepTape=F, temp=1, cache=new.env(hash =
    T, size = 10000L))}
\description{Undocumented.}
\details{method description}
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
\item{modejumping}{...}
\item{verbose}{...}
\item{keepTape}{...}
\item{temp}{...}
\item{cache}{...}}
