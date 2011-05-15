\name{as.grain}
\alias{as.grain}
\title{Convert to grain.}
\usage{as.grain(x, net, dat)}
\description{Convert to grain.}
\details{This is an interface to the package \pkg{gRain}. We supply a Bayesian
Network, its parameters, and the raw data. An object of class 'grain'
is returned. This object can be queried

Specifically, this function is a wrapper around
\code{\link[gRain]{cptable}}, \code{\link[gRain]{extractCPT}} and
\code{\link[gRain]{grain}}.}
\value{An object of class "grain".}
\seealso{\code{\link[gRain]{querygrain}} for querying the result. For
queries involving conditioning, use \code{\link[gRain]{setFinding}}
before running \code{querygrain}.}
\arguments{\item{x}{A list of multinomial parameters, in the form
output by \code{\link{ml}}}
\item{net}{a Bayesian Network. An object of class 'bn'}
\item{dat}{A data.frame with factors in columns. The columns should
correspond to the nodes of net, and the parameters of x.}}
