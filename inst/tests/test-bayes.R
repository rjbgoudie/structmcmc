# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Bayes parameters")

test_that("Bayes parameters (2)", {
  d <- data.frame(
    a = factor(c(1, rep(3, 2), rep(1, 7))),
    b = factor(c(2, rep(1, 2), 3, 3, rep(2, 5))),
    c = factor(c(2, rep(2, 3), rep(1, 6))),
    d = factor(c(1:3, 2:3, 1, 1, 3:2, 2))
  )

  net <- bn(integer(0), integer(0), integer(0), c(1, 2, 3))
  out <- bayes(net, d, prior = "qi")

  qi <- 1/(2 * 3 * 2)
  one <- (2 + qi)/(5 + 3 * qi)
  two <- (2 + qi)/(5 + 3 * qi)
  three <- (1 + qi)/(5 + 3 * qi)

  expect_that(out[[4]][["1,2,1"]], equals(c("1" = one,
                                            "2" = two,
                                            "3" = three)))
})
