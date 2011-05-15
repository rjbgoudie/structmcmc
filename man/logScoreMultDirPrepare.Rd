\name{logScoreMultDirPrepare}
\alias{logScoreMultDirPrepare}
\title{Internal functions.}
\usage{logScoreMultDirPrepare(data, logScoreParameters, checkInput=T)}
\description{Internal functions.}
\details{Convert a data frame to the appropriate format for the fast/incremental
logScoreMultDir functions, and return as part of the logScoreParameters
list.

In particular, the data is converted to a matrix, and the factor levels
taken integer values from 0, 1, .... ie not on 1, 2, 3.}
\value{A list with the contents of logScoreParameters, with the following
components added or altered:
data A matrix (NOT data.frame), with columns being integers
in the range 0, 1, 2, ....  giving the values of each
random variable.
nl   A numeric vector of length ncol(data), specifying the number of
levels that each random variable takes.}
\arguments{\item{data}{A data.frame, with columns being factors giving the
values of each random variable.}
\item{logScoreParameters}{A list, optionally containing other parameters.
In particular, it may contain 'hyperparameters'}
\item{checkInput}{A logical of length 1, specifying whether to check
the inputs to the function.}}
