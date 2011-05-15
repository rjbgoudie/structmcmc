\name{logScoreZellnerPrepare}
\alias{logScoreZellnerPrepare}
\title{Internal functions.}
\usage{logScoreZellnerPrepare(data, logScoreParameters, checkInput=T)}
\description{Internal functions.}
\details{Convert a data frame to the appropriate format for the fast/incremental
logScoreZellner functions, and return as part of the logScoreParameters
list.

In particular, the data is converted to a matrix, and the factor levels
taken integer values from 0, 1, .... ie not on 1, 2, 3.}
\value{A list with the contents of logScoreParameters, with the following
components added or altered:
\describe{
\item{data}{A matrix with columns giving the value of each random
variable.}
\item{nl}{A numeric vector of length ncol(data), specifying the number
of levels that each random variable takes.}
}}
\seealso{\code{\link{logScoreZellner}}}
\arguments{\item{data}{A data.frame, with columns being factors giving the values of
each random variable.}
\item{logScoreParameters}{A list with the following components:
\describe{
\item{data}{A matrix with columns giving the values of each random
variable.}
\item{nl}{A numeric vector of length nNodes(currentBN), specifying the
number of levels that each random variable takes.}
}}
\item{checkInput}{A logical of length 1, specifying whether to check the
inputs to the function.}}
