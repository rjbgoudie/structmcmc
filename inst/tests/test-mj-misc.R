# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("MJ Misc")

test_that("whichGraphsNotNeighbours", {
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  dat <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  p1 <- bn(integer(0), integer(0), 1)
  p2 <- bn(integer(0), 1, 1)
  p3 <- bn(integer(0), 1, integer(0))

  p1s <- logScoreMultDir(p1, dat)[[1]]
  p2s <- logScoreMultDir(p2, dat)[[1]]
  p3s <- logScoreMultDir(p3, dat)[[1]]

  # note these results have been checked for adequacy,
  # A better heuristic may come up with a better solution
  pl <- parental.list(p1, p1, p2)
  ms <- c(p1s, p1s, p2s)
  expect_that(whichGraphsNotNeighbours(pl, ms), equals(3))

  pl <- parental.list(p1, p1, p1, p1, p1, p2, p2)
  ms <- c(p1s, p1s, p1s, p1s, p1s, p2s, p2s)
  expect_that(whichGraphsNotNeighbours(pl, ms), equals(6))

  pl <- parental.list(p1, p2, p3)
  ms <- c(p1s, p2s, p3s)
  expect_that(whichGraphsNotNeighbours(pl, ms), equals(3))

  # flips
  p4 <- bn(2, integer(0), integer(0))
  p4s <- logScoreMultDir(p4, dat)[[1]]

  pl <- parental.list(p1, p3, p4)
  ms <- c(p1s, p3s, p4s)
  expect_that(whichGraphsNotNeighbours(pl, ms), equals(c(1, 3)))

  p5 <- bn(integer(0), 1, c(1, 2))
  p5s <- logScoreMultDir(p5, dat)[[1]]

  pl <- parental.list(p1, p2, p3, p4, p5)
  ms <- c(p1s, p2s, p3s, p4s, p5s)
  expect_that(whichGraphsNotNeighbours(pl, ms), equals(c(4, 5)))
})
