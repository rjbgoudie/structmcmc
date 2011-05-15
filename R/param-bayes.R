# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Bayesian posterior parameter estimates.
#'
#' A generic
#'
#' @param x An object
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{bayes.bn}}, \code{\link{ml}}
#' @examples
#' d <- data.frame(
#'   a = factor(c(1, rep(3, 2), rep(1, 7))),
#'   b = factor(c(2, rep(1, 2), 3, 3, rep(2, 5))),
#'   c = factor(c(2, rep(2, 3), rep(1, 6))),
#'   d = factor(c(1:3, 2:3, 1, 1, 3:2, 2))
#' )
#' 
#' net <- bn(integer(0), integer(0), integer(0), c(1, 2, 3))
#' bayes(net, d, prior = "qi")
bayes <- function (x, ...) {
  UseMethod("bayes")
}

#' Bayesian posterior parameter estimates.
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
#' @examples
#' d <- data.frame(
#'   a = factor(c(1, rep(3, 2), rep(1, 7))),
#'   b = factor(c(2, rep(1, 2), 3, 3, rep(2, 5))),
#'   c = factor(c(2, rep(2, 3), rep(1, 6))),
#'   d = factor(c(1:3, 2:3, 1, 1, 3:2, 2))
#' )
#' 
#' net <- bn(integer(0), integer(0), integer(0), c(1, 2, 3))
#' bayes(net, d, prior = "qi")
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
