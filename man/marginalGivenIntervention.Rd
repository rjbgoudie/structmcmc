\name{marginalGivenIntervention}
\alias{marginalGivenIntervention}
\title{Get marginal probabilities, given an intervention.}
\usage{marginalGivenIntervention(net, node, dat, FUN=identity, ...)}
\description{Get marginal probabilities, given an intervention.}
\details{Get the marginal probability of a node, when each other node (separately)
has been conditioned upon each of its levels.}
\value{A list, with a component for each node of the Bayes Net (apart
from the node \code{node}). Each component is a list, corresponding to
the named level. These components contain the result of applying
result of applying \code{FUN} to the result of \code{querygrain}.}
\arguments{\item{net}{A BN. A "parental" object.}
\item{node}{An integer or a name of a column of \code{dat}}
\item{dat}{A data frame, with columns corresponding to the Bayes Net
in \code{grain}.}
\item{FUN}{a function that is applied to the output of
\code{\link[gRain]{querygrain}}. The first argument should accept
the result of \code{\link[gRain]{querygrain}} FOR A PARTICULAR NODE.
(Usaully \code{querygrain} returns a list, since it may provide the
distribution of a number of nodes at once. Here, we remove that list.)
Any other arguments are passed \code{...}}
\item{...}{Passed to \code{FUN}}}
