\name{routesRemoveEdges}
\alias{routesRemoveEdges}
\title{Update a routes matrix (edge removal)...}
\usage{routesRemoveEdges(x, i, j)}
\description{Update a routes matrix (edge removal)}
\details{A routes matrix is a matrix A, such that each element (i, j) is the
number of routes from i to j in some directed graph.

This function updates the routes matrix to account for the deletion of
an edge from i to j in the directed graph}
\arguments{\item{x}{A routes matrix}
\item{i}{The nodes from which the removed edge emanates}
\item{j}{The node that the removed edge goes to. ONLY ONE NODE HERE}}
