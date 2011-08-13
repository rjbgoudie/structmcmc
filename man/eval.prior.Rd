\name{eval.prior}
\alias{eval.prior}
\title{Compute prior score of a Bayesian network.}
\usage{eval.prior(x, localPriors)}

\description{
  Compute prior score of a Bayesian network.
}

\details{
  Computes the complete prior of a network, from a
  \code{localPriors} object
}
\value{The prior of the complete network \code{x}}
\arguments{
  \item{x}{An object of class 'bn'.}
  \item{localPriors}{A list of functions of the same length as \code{x}
that returns the local prior score of the corresponding node, given a
numeric vector of parents.}
}
