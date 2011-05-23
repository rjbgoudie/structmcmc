\name{mukherjeeBioinformaticsPrior}
\alias{mukherjeeBioinformaticsPrior}
\title{Mukherjee Bioinformatics prior.}
\usage{mukherjeeBioinformaticsPrior(x, k0, kmax, lambda)}

\description{
  Mukherjee Bioinformatics prior.
}

\details{
  Compute prior score P(M) for a 'bvsresponse'. The prior
  is flat from 0 parents to k0 parents, then
  exponentially(lambda) decreaasing to kmax, at which point
  it becomes zero.

  The prior originates in Mukherjee el al. (2009), and is
  equation (8) in section 3.1.3.

  Mukherjee et al. Sparse combinatorial inference with an
  application in cancer biology. Bioinformatics (2009) vol.
  25 (2) pp. 265
}
\value{The prior score for the 'bvsresponse' x}
\arguments{
  \item{x}{A 'bvsresponse' object}
  \item{k0}{The level below which the prior stops being flat}
  \item{kmax}{The level above which the prior is zero}
  \item{lambda}{The exponential decay parameter}
}
