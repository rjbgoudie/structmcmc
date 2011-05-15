\name{parentalFromEPThreshold}
\alias{parentalFromEPThreshold}
\title{Threshold posterior edge probabilities.}
\usage{parentalFromEPThreshold(ep, threshold)}
\description{Threshold posterior edge probabilities.}
\details{Return the parental given by thresholding a edge probability matrix
at a given level. The inequality is >=. Note this may well be cyclic.}
\value{An object of class parental}
\arguments{\item{ep}{A matrix of class ep}
\item{threshold}{The value at which to threshold the matrix}}
