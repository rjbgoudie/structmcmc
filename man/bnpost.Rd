\name{bnpost}
\alias{bnpost}
\title{method name...}
\usage{bnpost(bnspace, logScore, data, logScoreFUN=logScoreMultDir)}
\description{method name}
\details{method description}
\arguments{\item{bnspace}{...}
\item{logScore}{...}
\item{data}{...}
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
}}}
