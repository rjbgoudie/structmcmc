# Part of the "structural" package, http://github.com/rbtgde/structural
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structural
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Bayesian posterior parameter estimates
#'
#' A generic
#'
#' @param x An object
#' @param ... Further arguments passed to method
#' @export
bayes <- function (x, ...) {
  UseMethod("bayes")
}

#' Bayesian posterior parameter estimates
#'
#' A wrapper around \code{\link{ml}}.
#'
#' This is the probability. This is also the expectation. See
#' Neapolitean p379.
#'
#' @param x The Bayesian Network. An object of class 'bn'
#' @param data A data frame
#' @param nodes A subset of 1, ..., \code{nNodes(x)}. A numeric vector.
#' @param prior Only \code{"qi"} is implemented at the moment
#' @param ... Further arguments, passed to \code{\link{ml}}
#' @return As \code{\link{ml}}
#' @S3method bayes bn
#' @method bayes bn
bayes.bn <- function(x,
                  data,
                  nodes = seq_along(x),
                  prior = "qi",

                  ...){
  ml(x              = x,
     data           = data,
     nodes          = nodes,
     regularisation = prior,
     ...)
}