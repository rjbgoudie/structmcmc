# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Cross-validation tests")

test_that("3-node Bayesian Network", {
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  train <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  x <- bn(integer(0), c(1,3), integer(0))

  test <- data.frame(x1 = x1[c(2,3,4,2,3)],
                     x2 = x2[c(6,4,2,4,2)],
                     x3 = x3[c(2,3,4,3,4)])

  bayes(x, train)

  residualsMultDir(x, train, test)
})

test_that("3-node Bayesian Network", {
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  train <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  x <- bn(integer(0), c(3), integer(0))
  y <- bn(integer(0), c(1), integer(0))

  test <- data.frame(x1 = x1[c(2,3,4,2,3,8,7)],
                     x2 = x2[c(6,4,2,4,2,7,8)],
                     x3 = x3[c(2,3,4,3,4,7,7)])

  predictModelAverageNode(2, list(x, y), train, test, weights = c(0.1, 0.9))

  expect_that(residualsMultDir(bn.list(x, y), weights = c(0.1, 0.9), train, test), equals(c(4, 4, 0))
})
