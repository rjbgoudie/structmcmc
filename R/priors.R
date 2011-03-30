# Part of the "structural" package, http://github.com/rjbgoudie/structural
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rjbgoudie/structural
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Compute prior score P(M) for a 'bvsresponse'. The prior is flat from
#' 0 parents to k0 parents, then exponentially(lambda) decreaasing to kmax,
#' at which point it becomes zero.
#'
#' The prior originates in Mukherjee el al. (2009), and is equation (8) in
#' section 3.1.3.
#'
#' Mukherjee et al. Sparse combinatorial inference with an application in
#' cancer biology. Bioinformatics (2009) vol. 25 (2) pp. 265
#'
#' @param x A 'bvsresponse' object
#' @param k0 The level below which the prior stops being flat
#' @param kmax The level above which the prior is zero
#' @param lambda The exponential decay parameter
#' @return The prior score for the 'bvsresponse' x
#' @export
mukherjeeBioinformaticsPrior <- function(x, k0, kmax, lambda){
  stopifnot("bvsresponse" %in% class(x))

  numberOfParents <- length(x$parents)
  kmaxTest <- numberOfParents <= kmax
  if (kmaxTest){
    part <- min(0, k0 - numberOfParents)
    exp(lambda * part)
  }
  else {
    0
  }
}

#' Checks the output of a LOG prior for validity. Basically is it positive?
#'
#' @param x A numeric of length 1. The LOG output of a prior function.
#' @return A logical of length 1, indicating the validity of x
#' @S3method is.valid prior
#' @method is.valid prior
#' @export
is.valid.prior <- function(x){
  tryCatch({
    all(class(x) %in% c("numeric", "integer"),
        !is.nan(x),
        is.finite(x))
    },
    error = function(e) F)
}
