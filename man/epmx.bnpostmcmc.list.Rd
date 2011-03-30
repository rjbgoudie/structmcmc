\name{epmx.bnpostmcmc.list}
\alias{epmx.bnpostmcmc.list}
\title{Computes the edge probabilities and return a matrix with these.}
\description{Computes the edge probabilities and return a matrix with these. The
format of the matrix is designed for the plotting function xyplot.epmx.}
\details{For a problem with k nodes, the output will have k^2 columns and nbin
rows. Columns are in order 1->1, 1->2, 1->3, ...., 2->1, 2->2 etc}
\value{An object of class "epmx", a matrix of the form described above.}
\arguments{\item{x}{An object of class "bnpostmcmc.list"}
\item{nbin}{The number of equally-sized bins into which the samples are
divided into before computing the edge probabilities of each}
\item{start}{An integer of length 1, specifying the amount of burn-in.
The samples start:end inclusive will be used.}
\item{end}{An integer of length 1, specifying the number of samples at the
end to ignore. The samples start:end inclusive will be used.}
\item{...}{Further arguments (unused)}}
