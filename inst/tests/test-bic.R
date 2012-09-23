# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Bayesian Information Criterion")

test_that("1 Node input", {
  d <- data.frame(
    a = factor(c(rep(1, 5), rep(2, 3), rep(1, 2), 2)),
    b = factor(c(rep(1, 5), rep(1, 3), rep(2, 2), 2))
  )

  net <- bn(integer(0), c(1))
  out <- scoreBIC(net, d)

  hand <- log((7/11)^7 * (4/11)^4 * (5/7)^5 * (2/7)^2 * (3/4)^3 * (1/4)^1) -
          3/2*log(11)

  expect_that(out, equals(hand))
})
