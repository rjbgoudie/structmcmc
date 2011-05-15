# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Setup constraint.
#'
#' ...
#'
#' @param constraint ..
#' @param initial ..
#' @export
#' @seealso \code{\link{satisfiesConstraint}},
#'   \code{\link{enforceForbiddenConstraint}}
setupConstraint <- function(constraint, initial){
  numberOfNodes <- nNodes(initial)
  if (is.null(constraint)){
    matrix(0, numberOfNodes, numberOfNodes)
  }
  else {
    stopifnot(class(constraint)    ==   "matrix",
              all(constraint       %in% c(-1, 0, 1)),
              all(diag(constraint) ==   0))

    # check initial meet constraint
    if (!satisfiesConstraint(initial, constraint)){
      stop("Initial network does not satisfy constraint")
    }
    constraint
  }
}

#' Undocumented.
#'
#' method description
#'
#' @param x ...
#' @param constraint ...
#' @export
#' @seealso \code{\link{setupConstraint}},
#'   \code{\link{enforceForbiddenConstraint}}
satisfiesConstraint <- function(x, constraint){
  stopifnot("parental" %in% class(x),
            class(constraint) == "matrix",
            all(constraint %in% c(-1, 0, 1)),
            all(diag(constraint) == 0))

  edgesRequiredByConstraint <- which(constraint == 1, arr.ind = T)
  edgesForbiddenByConstraint <- which(constraint == -1, arr.ind = T)
  satisfiesConstraint <- T

  for (row in seq_len(nrow(edgesRequiredByConstraint))){
    i <- edgesRequiredByConstraint[row, 1]
    j <- edgesRequiredByConstraint[row, 2]
    if (!i %in% x[[j]]){
      satisfiesConstraint <- F
    }

    # forbid the reverse of required link
    constraint[j, i] <- -1
  }

  for (row in seq_len(nrow(edgesForbiddenByConstraint))){
    i <- edgesForbiddenByConstraint[row, 1]
    j <- edgesForbiddenByConstraint[row, 2]
    if (i %in% x[[j]]){
      satisfiesConstraint <- F
    }
  }
  satisfiesConstraint
}

#' Undocumented.
#'
#' method description
#'
#' @param x ...
#' @param constraint ...
#' @export
#' @seealso \code{\link{satisfiesConstraint}}, \code{\link{setupConstraint}}
enforceForbiddenConstraint <- function(x, constraint){
  stopifnot("parental" %in% class(x),
            class(constraint) == "matrix",
            all(constraint %in% c(-1, 0, 1)),
            all(diag(constraint) == 0))

  edgesForbiddenByConstraint <- which(constraint == -1, arr.ind = T)

  for (row in seq_len(nrow(edgesForbiddenByConstraint))){
    i <- edgesForbiddenByConstraint[row, 1]
    j <- edgesForbiddenByConstraint[row, 2]
    x[[j]] <- setdiff(x[[j]], i)
  }
  x
}


#' Resample a pair of nodes together.
#'
#' 
#'
#' @param constraint ..
#' @export
getRequiredFromConstraint <- function(constraint){
  numberOfNodes <- dim(constraint)[1]
  nodesSeq <- seq_len(numberOfNodes)
  lapply(nodesSeq, whichNum, x = constraint, i = 1)
}

#' Resample a pair of nodes together.
#'
#' 
#'
#' @param constraint ..
#' @export
getBannedFromConstraint <- function(constraint){
  numberOfNodes <- dim(constraint)[1]
  nodesSeq <- seq_len(numberOfNodes)
  lapply(nodesSeq, whichNum, x = constraint, i = -1)
}
