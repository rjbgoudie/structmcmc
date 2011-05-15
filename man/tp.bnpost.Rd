\name{tp.bnpost}
\alias{tp.bnpost}
\title{Transition probabilities.}
\usage{\method{tp}{bnpost} (x, sampler="mh", allowFlips=T, verbose=F, ...)}
\description{Transition probabilities.}
\details{Computes the matrix transition probabilities for the specified sampler.}
\value{A matrix of transition probabilities.}
\seealso{\code{\link{tp}}}
\arguments{\item{x}{An object of class "bnpost"}
\item{sampler}{Which sampler to use. Only "mh" for Metropolis-Hastings
is implemented}
\item{allowFlips}{A logical of length 1, specifying whether the sampler is
allowed to reverse the direction of single edges?}
\item{verbose}{...}
\item{...}{Further arguments (unused)}}
