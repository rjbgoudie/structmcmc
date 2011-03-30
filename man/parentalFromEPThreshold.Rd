\name{parentalFromEPThreshold}
\alias{parentalFromEPThreshold}
\title{Return the parental given by thresholding a edge probability matrix...}
\usage{parentalFromEPThreshold(ep, threshold)}
\description{Return the parental given by thresholding a edge probability matrix
at a given level. The inequality is >=. Note this may well be cyclic.}
\value{An object of class parental}
\arguments{\item{ep}{A matrix of class ep}
\item{threshold}{The value at which to threshold the matrix}}
