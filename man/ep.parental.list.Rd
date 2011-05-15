\name{ep.parental.list}
\alias{ep.parental.list}
\title{Posterior edge probabilities.}
\usage{\method{ep}{parental.list} (x, nbin=1, start=1, end=length(x),
    verbose=F, ...)}
\description{Posterior edge probabilities.}
\details{Calculate edge probabilities given a list of
parental graphs (of class parental.list)}
\value{if nbin == 1:
A matrix of class 'ep' with entry (i,j) containing the probability of
an edge from node i --> j
if nbin > 1:
A list of class ep.list, containing matrices as described above for
each of the nbin bins into which the parental.list was split}
\arguments{\item{x}{A parental.list}
\item{nbin}{The number of equally-sized bins into which parental.list into
before computing the edge probabilities of each}
\item{start}{An integer of length 1, specifying the amount of burn-in.
The samples start:end inclusive will be used.}
\item{end}{An integer of length 1, specifying the number of samples at the
end to ignore. The samples start:end inclusive will be used.}
\item{verbose}{...}
\item{...}{Further arguments (unused)}}
