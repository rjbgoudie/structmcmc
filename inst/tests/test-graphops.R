# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Graph ops")

test_that("getPossibleParents", {
  bn <- bn(3:4, c(), c(), c())
  nonDescendantsList <- list(2:4, c(1L, 3L, 4L), c(1L, 4L), 1:3)
  descendantsList <- list(1L, 2L, 2:3, 4L)
  numberOfNodes <- 4
  change <- c(1, 3, 4)
  maxIndegree <- 4

  actual <- getPossibleParents(bn                 = bn,
                               nonDescendantsList = nonDescendantsList,
                               descendantsList    = descendantsList,
                               numberOfNodes      = numberOfNodes,
                               change             = change,
                               maxIndegree        = maxIndegree)
  expected <- list(c(2L, 3L, 4L), NULL, integer(0), integer(0))
  expect_identical(actual, expected)

  bn <- bn(c(), c(), c(1L, 4L), c())
  actual <- getPossibleParents(bn                 = bn,
                               nonDescendantsList = nonDescendantsList,
                               descendantsList    = descendantsList,
                               numberOfNodes      = numberOfNodes,
                               change             = change,
                               maxIndegree        = maxIndegree)
  expected <- list(integer(0), NULL, c(1L, 4L), integer(0))
  expect_identical(actual, expected)
})
