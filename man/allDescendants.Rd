\name{allDescendants}
\alias{allDescendants}
\title{Get descendants of all nodes.}
\usage{allDescendants(currentNetwork)}
\description{Get descendants of all nodes.}
\details{Note that the descendants of each node includes that node!}
\value{A list, each componentof which contains a list of the
descendants of the corresponding node of the supplied BN.}
\arguments{\item{currentNetwork}{A list, containing in the first position the
starting \code{bn}, and in the second position the routes matrix for
that BN.}}
