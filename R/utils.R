# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Check validity of start and end.
#'
#' Checks if start and end are valid as start and end points
#' for taking a window (a subset, e.g. to remove burn-in) in a MCMC run
#'
#' @param start The start value to check
#' @param end The end value to check
#' @param length A numeric of length 1. The length of the MCMC run that start
#'   and end should window.
#' @return A logical of length 1 indicating if start and end are valid
validStartEnd <- function(start, end, length){
  all(isTRUE(is.wholenumber(start)),
      isTRUE(is.wholenumber(end)),
      end <= length,
      start < end,
      start >= 1
  )
}

#' Fast, dangerous row sums.
#'
#' A fast, simple version of \code{rowSums}.
#' This version only handles matrices.
#'
#' Note that no sanity checks are performed on the input.
#'
#' @param x A matrix.
#' @return A vector of row sums.
rowSums2 <- function(x){
  dn <- dim(x)
  p <- prod(dn[-(1L)])
  dn <- dn[1L]
  z <- .Internal(rowSums(x, prod(dn), p, na.rm = F))
  if (length(dn) > 1L) {
      dim(z) <- dn
  }
  z
}

#' Fast, dangerous set difference.
#'
#' A fast, simple version of \code{setdiff}.
#' This version does not handle factors.
#'
#' Note that no sanity checks are performed on the input.
#'
#' @param x A vector, of the same mode as \code{y}.
#' @param y A vector, of the same mode as \code{x}.
#' @return A vector of the same mode as the inputs.
setdiff3 <- function(x, y){
  x[fastmatch::fmatch(x, y, 0L, NULL) == 0L]
}


#' Fast, dangerous set intersect.
#'
#' A fast, simple version of \code{intersect}.
#' This version does not handle factors.
#'
#' Note that no sanity checks are performed on the input.
#'
#' @param x A vector, of the same mode as \code{y}.
#' @param y A vector, of the same mode as \code{x}.
#' @return A vector of the same mode as the inputs.
intersect2 <- function(x, y, nmax = NA){
  .Internal(unique(x             = y[fastmatch::fmatch(x, y, 0L, NULL)],
                   incomparables = F,
                   fromLast      = F,
                   nmax          = nmax))
}

#' Fast, dangerous matrix generation.
#'
#' A fast, simple version of \code{matrix}.
#'
#' The matrix must be supplied by column (i.e. \code{byrow = FALSE}) and
#' the matrix can not have any names (i.e. \code{dimnames = NULL}).
#'
#' Note that no sanity checks are performed on the input.
#'
#' @param data A data vector.
#' @param nrow The desired number of rows
#' @param ncol The desired number of columns
#' @return A matrix of dimension \code{nrow} by \code{ncol}, containing
#'   data \code{data}.
matrix2 <- function(data, nrow, ncol)
{
  .Internal(matrix(data, nrow, ncol, FALSE, NULL,
                   missing(nrow), missing(ncol)))
}

#' Find rows satisfying equality.
#'
#' Find which rows in column \code{col} of a matrix \code{x} are equal to
#' \code{i}.
#'
#' Note that no sanity checks are performed on the input.
#'
#' @param col The column of the matrix to use. A numeric vector of length 1.
#' @param x A matrix.
#' @param i The value to check for
#' @return A numeric vector, of the
#' @export
whichNum <- function(col, x, i){
  which(x[, col] == i)
}
