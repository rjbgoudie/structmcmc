# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Priors")

test_that("Sach's Bioinformatics Prior", {
  net <- as.bvsresponse(bvs(integer(0), integer(0)), response = 2)
  test1 <- mukherjeeBioinformaticsPrior(x      = net,
                                        k0     = 1,
                                        kmax   = 1,
                                        lambda = 2)

  expect_that(test1, equals(1))

  net <- as.bvsresponse(bvs(integer(0), 1), response = 2)
  test2 <- mukherjeeBioinformaticsPrior(x      = net,
                                        k0     = 1,
                                        kmax   = 1,
                                        lambda = 2)

  expect_that(test2, equals(1))

  net <- as.bvsresponse(bvs(integer(0), integer(0), 1L))
  test3 <- mukherjeeBioinformaticsPrior(x      = net,
                                        k0     = 1,
                                        kmax   = 2,
                                        lambda = 2)

  expect_that(test3, equals(1))

  net <- as.bvsresponse(bvs(integer(0), integer(0), c(1L, 2L)))
  test4 <- mukherjeeBioinformaticsPrior(x      = net,
                                        k0     = 1,
                                        kmax   = 2,
                                        lambda = 2)

  expect_that(test4, equals(exp(-2)))

  net <- as.bvsresponse(bvs(integer(0), integer(0), c(2L, 4L, 5L, 6L, 7L),
             integer(0), integer(0), integer(0),
             integer(0), integer(0), integer(0)))
  test5 <- mukherjeeBioinformaticsPrior(x      = net,
                                        k0     = 3,
                                        kmax   = 7,
                                        lambda = 0.5)

  expect_that(test5, equals(exp(0.5 * -2)))
})