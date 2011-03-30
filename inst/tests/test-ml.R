# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("ML Parameters")

test_that("Very weird input", {

})

test_that("1 Node input", {
  d <- data.frame(
    a = factor(c(1,2))
  )

  net <- bn(integer(0))
  out <- ml(net, d)
  expect_that(out[[1]][[1]], equals(c("1" = 0.5, "2" = 0.5)))
})

test_that("3 Node input", {
  d <- data.frame(
    a = factor(c(1, rep(3,2), rep(1, 7))),
    b = factor(c(2, rep(1, 4), rep(2, 5))),
    c = factor(c(2, rep(2, 3), rep(1, 6)))
  )

  net <- bn(integer(0), integer(0), c(1,2))
  out <- ml(net, d)

  expect_that(out[[1]][[1]], equals(c("1" = 0.8, "3" = 0.2)))
  expect_that(out[[2]][[1]], equals(c("1" = 0.4, "2" = 0.6)))

  expect_that(out[[3]][["1,1"]], equals(c("1" = 0.5, "2" = 0.5)))
  expect_that(out[[3]][["1,2"]], equals(c("1" = 5/6, "2" = 1/6)))
  expect_that(out[[3]][["3,1"]], equals(c("1" = 0, "2" = 1)))
  expect_that(out[[3]][["3,2"]], equals(c("1" = NaN, "2" = NaN)))
})

test_that("NaNs", {
  d <- data.frame(
    a = factor(c(1, rep(3,2), rep(1, 7))),
    b = factor(c(2, rep(1, 4), rep(2, 5))),
    c = factor(c(2, rep(2, 3), rep(1, 6)))
  )

  net <- bn(integer(0), integer(0), c(1,2))
  out <- ml(net, d, regularisation = 0)

  expect_that(out[[3]][["3,2"]], equals(c("1" = NaN, "2" = NaN)))
})

test_that("qi regularisation", {
  d <- data.frame(
    a = factor(c(1, rep(3,2), rep(1, 7))),
    b = factor(c(2, rep(1, 4), rep(2, 5))),
    c = factor(c(2, rep(2, 3), rep(1, 6)))
  )

  net <- bn(integer(0), integer(0), c(1,2))
  out <- ml(net, d, regularisation = "qi")

  expect_that(out[[3]][["1,1"]], equals(c("1" = 0.5, "2" = 0.5)))
  expect_that(out[[3]][["1,2"]], equals(c("1" = (5 + 1/(2 * 2))/(6 + 2 * 1/(2 * 2)),
                                          "2" = (1 + 1/(2 * 2))/(6 + 2 * 1/(2 * 2)))))
  expect_that(out[[3]][["3,1"]], equals(c("1" = (0 + 1/(2 * 2))/(2 + 2 * 1/(2 * 2)),
                                          "2" = (2 + 1/(2 * 2))/(2 + 2 * 1/(2 * 2)))))
  expect_that(out[[3]][["3,2"]], equals(c("1" = 0.5, "2" = 0.5)))
})

test_that("regularisation = 1", {
  d <- data.frame(
    a = factor(c(1, rep(3,2), rep(1, 7))),
    b = factor(c(2, rep(1, 4), rep(2, 5))),
    c = factor(c(2, rep(2, 3), rep(1, 6)))
  )

  net <- bn(integer(0), integer(0), c(1,2))
  out <- ml(net, d, regularisation = 1)

  expect_that(out[[3]][["1,1"]], equals(c("1" = 0.5, "2" = 0.5)))
  expect_that(out[[3]][["1,2"]], equals(c("1" = (5 + 1)/(6 + 2 * 1),
                                          "2" = (1 + 1)/(6 + 2 * 1))))
  expect_that(out[[3]][["3,1"]], equals(c("1" = (0 + 1)/(2 + 2 * 1),
                                          "2" = (2 + 1)/(2 + 2 * 1))))
  expect_that(out[[3]][["3,2"]], equals(c("1" = 0.5, "2" = 0.5)))
})

test_that("regularisation = 3", {
  d <- data.frame(
    a = factor(c(1, rep(3,2), rep(1, 7))),
    b = factor(c(2, rep(1, 4), rep(2, 5))),
    c = factor(c(2, rep(2, 3), rep(1, 6)))
  )

  net <- bn(integer(0), integer(0), c(1,2))
  out <- ml(net, d, regularisation = 3)

  expect_that(out[[3]][["1,1"]], equals(c("1" = 0.5, "2" = 0.5)))
  expect_that(out[[3]][["1,2"]], equals(c("1" = (5 + 3)/(6 + 2 * 3),
                                          "2" = (1 + 3)/(6 + 2 * 3))))
  expect_that(out[[3]][["3,1"]], equals(c("1" = (0 + 3)/(2 + 2 * 3),
                                          "2" = (2 + 3)/(2 + 2 * 3))))
  expect_that(out[[3]][["3,2"]], equals(c("1" = 0.5, "2" = 0.5)))
})

test_that("qi regularisation (2)", {
  d <- data.frame(
    a = factor(c(1, rep(3,2), rep(1, 7))),
    b = factor(c(2, rep(1, 2), 3, 3, rep(2, 5))),
    c = factor(c(2, rep(2, 3), rep(1, 6))),
    d = factor(c(1:3, 2:3, 1, 1, 3:2, 2))
  )

  net <- bn(integer(0), integer(0), integer(0), c(1, 2, 3))
  out <- ml(net, d, regularisation = "qi")

  qi <- 1/(2 * 3 * 2)
  one <- (2 + qi)/(5 + 3 * qi)
  two <- (2 + qi)/(5 + 3 * qi)
  three <- (1 + qi)/(5 + 3 * qi)

  expect_that(out[[4]][["1,2,1"]], equals(c("1" = one,
                                            "2" = two,
                                            "3" = three)))
})