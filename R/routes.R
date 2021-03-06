# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Update a routes matrix (edge addition).
#'
#' A routes matrix is a matrix A, such that each element (i, j) is the
#' number of routes from i to j in some directed graph.
#'
#' This function updates the routes matrix to account for the addition of
#' an edge from i to j in the directed graph
#'
#' @param x A routes matrix
#' @param i The nodes from which the added edge emanates
#' @param j The node that the added edge goes to. ONLY ONE NODE HERE
#' @export
#' @seealso \code{\link{routesRemoveEdges}}, \code{\link{nonDescendants}}
routesAddEdges <- function(x, i, j){
  if (length(i) > 0){
    x + rowSums2(x[, i, drop = F]) * x[rep(j, dim(x)[1]), ]
  } else {
    x
  }
}

#' Update a routes matrix (edge removal).
#'
#' A routes matrix is a matrix A, such that each element (i, j) is the
#' number of routes from i to j in some directed graph.
#'
#' This function updates the routes matrix to account for the deletion of
#' an edge from i to j in the directed graph
#'
#' @param x A routes matrix
#' @param i The nodes from which the removed edge emanates
#' @param j The node that the removed edge goes to. ONLY ONE NODE HERE
#' @export
#' @seealso \code{\link{routesAddEdges}}, \code{\link{nonDescendants}}
routesRemoveEdges <- function(x, i, j){
  if (length(i) > 0){
    x - rowSums2(x[, i, drop = F]) * x[rep(j, dim(x)[1]), ]
  } else {
    x
  }
}

#' Find nonDescendants.
#'
#' Find the non-descendants of a node
#'
#' @param x A routes matrix
#' @param node A node
#' @export
#' @seealso \code{\link{routesAddEdges}}, \code{\link{routesRemoveEdges}}
nonDescendants <- function(x, node){
  .Internal(which(x[node, ] == 0))
}
