# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

library(testthat)
library(structmcmc)

if (as.numeric(R.version$minor) < 11){
  vapply <- function(X, FUN, FUN.VALUE, ..., USE.NAMES = TRUE){
  sapply(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE)
  }
}

is_within <- function(expected, tolerance){
  name <- deparse(substitute(expected))
  function(actual) {
      testthat:::expectation(
        identical(all.equal.numeric(expected, actual, tolerance = tolerance, scale = 1), TRUE),
        paste("was out by ", abs(expected - actual), " when max is ", tolerance, sep = "")
      )
    }
}

test_package("structmcmc")
