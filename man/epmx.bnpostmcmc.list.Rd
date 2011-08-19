\name{epmx.bnpostmcmc.list}
\alias{epmx.bnpostmcmc.list}
\title{Edge probabilities matrix.}
\usage{
  \method{epmx}{bnpostmcmc.list} (x, method = "online",
  nbin = floor((end - start + 1)/100), start = 1, end =
  length(x[[1]]), verbose = F, ...)
}
\arguments{
  \item{x}{An object of class "bnpostmcmc.list"}

  \item{method}{Either \code{"online"}, or
  \code{"offline"}. For online, the edge totals kept online
  will be used. This}

  \item{nbin}{The number of equally-sized bins into which
  the samples are divided into before computing the edge
  probabilities of each. This applies only for \code{method
  = "offline"}.}

  \item{start}{An integer of length 1, specifying the
  amount of burn-in. The samples start:end inclusive will
  be used.}

  \item{end}{An integer of length 1, specifying the number
  of samples at the end to ignore. The samples start:end
  inclusive will be used.}

  \item{verbose}{Should progress be shown? A logical.}

  \item{...}{Further arguments (unused)}
}
\value{
  An object of class "epmx", a matrix of the form described
  above.
}
\description{
  Computes the edge probabilities and return a matrix with
  these. The format of the matrix is designed for the
  plotting function splom.epmx.
}
\details{
  For a problem with k nodes, the output will have k^2
  columns and nbin rows. Columns are in order 1->1, 1->2,
  1->3, ...., 2->1, 2->2 etc
}
\seealso{
  \code{\link{epmx}}, \code{\link{splom.epmx}},
  \code{\link{splom.bnpostmcmc.list}}
}

