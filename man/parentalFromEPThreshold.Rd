\name{parentalFromEPThreshold}
\alias{parentalFromEPThreshold}
\title{Threshold posterior edge probabilities.}
\usage{parentalFromEPThreshold(ep, threshold)}
\description{Threshold posterior edge probabilities.}
\details{Return the parental given by thresholding a edge probability matrix
at a given level. The inequality is >=. Note this may well be cyclic.}
\value{An object of class parental}
\seealso{\code{\link{ep}}}
\arguments{\item{ep}{A matrix of class ep}
\item{threshold}{The value at which to threshold the matrix}}
\examples{x1 <- factor(c(1, 1, 0, 1))
x2 <- factor(c(0, 1, 0, 1))
x <- data.frame(x1 = x1, x2 = x2)

exact <- posterior(data = x, "exact")

myep <- ep(exact)
parentalFromEPThreshold(myep, 0.2)
parentalFromEPThreshold(myep, 0.4)}
