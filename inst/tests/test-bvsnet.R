# Part of the "structural" package, http://github.com/rjbgoudie/structural
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rjbgoudie/structural
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

# context("BVS net")
#
# test_that("2-node Bayesian Network", {
#   x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
#   x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
#   x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
#   dat <- data.frame(x1 = x1, x2 = x2,  x3 = x3)
#  
#   set.seed(10)
#   out <- bvsnet(dat,
#                 nIterations = 10000,
#                 nBurnIn     = 2500,
#                 threshold   = 0.7,
#                 verbose     = F)
#  
#   expected <-
#   structure(list(structure(list(integer(0), 3L, integer(0)), class = c("bn",
#   "parental")), structure(list(integer(0), 3L, integer(0)), class = c("bn",
#   "parental")), structure(list(integer(0), 3L, integer(0)), class = c("bn",
#   "parental")), structure(list(integer(0), 3L, integer(0)), class = c("bn",
#   "parental")), structure(list(integer(0), integer(0), 2L), class = c("bn",
#   "parental")), structure(list(integer(0), integer(0), 2L), class = c("bn",
#   "parental")), structure(list(integer(0), integer(0), 2L), class = c("bn",
#   "parental")), structure(list(integer(0), integer(0), 2L), class = c("bn",
#   "parental")), structure(list(integer(0), 3L, integer(0)), class = c("bn",
#   "parental")), structure(list(integer(0), 3L, integer(0)), class = c("bn",
#   "parental"))), class = c("bn.list", "parental.list"))
#  
#   # expected <- list(structure(list(integer(0), 3L, integer(0)),
#   #                            class = c("bn", "parental")),
#   #                  structure(list(integer(0), integer(0), 2L),
#   #                            class = c("bn", "parental")))
#  
#   expect_that(out, is_identical_to(expected))
#  
#   set.seed(10)
#   out <- bvsnet(dat,
#                 nOrders     = 20,
#                 nIterations = 10000,
#                 nBurnIn     = 2500,
#                 threshold   = 0.1,
#                 verbose     = F)
#  
#   expected <-
#   structure(list(structure(list(integer(0), c(1L, 3L), 1L), class = c("bn",
#   "parental")), structure(list(integer(0), c(1L, 3L), 1L), class = c("bn",
#   "parental")), structure(list(3L, c(1L, 3L), integer(0)), class = c("bn",
#   "parental")), structure(list(2:3, 3L, integer(0)), class = c("bn",
#   "parental")), structure(list(2L, integer(0), 1:2), class = c("bn",
#   "parental")), structure(list(2L, integer(0), 1:2), class = c("bn",
#   "parental")), structure(list(2:3, integer(0), 2L), class = c("bn",
#   "parental")), structure(list(2L, integer(0), 1:2), class = c("bn",
#   "parental")), structure(list(3L, c(1L, 3L), integer(0)), class = c("bn",
#   "parental")), structure(list(integer(0), c(1L, 3L), 1L), class = c("bn",
#   "parental")), structure(list(3L, c(1L, 3L), integer(0)), class = c("bn",
#   "parental")), structure(list(2:3, 3L, integer(0)), class = c("bn",
#   "parental")), structure(list(integer(0), c(1L, 3L), 1L), class = c("bn",
#   "parental")), structure(list(3L, c(1L, 3L), integer(0)), class = c("bn",
#   "parental")), structure(list(3L, c(1L, 3L), integer(0)), class = c("bn",
#   "parental")), structure(list(integer(0), c(1L, 3L), 1L), class = c("bn",
#   "parental")), structure(list(2L, integer(0), 1:2), class = c("bn",
#   "parental")), structure(list(integer(0), c(1L, 3L), 1L), class = c("bn",
#   "parental")), structure(list(integer(0), c(1L, 3L), 1L), class = c("bn",
#   "parental")), structure(list(integer(0), 1L, 1:2), class = c("bn",
#   "parental"))), class = c("bn.list", "parental.list"))
#  
#   # this used to be the expectation. it seems wrong to me now?
#   # expected <- list(structure(list(integer(0), c(1L, 3L), 1L),
#   #                            class = c("bn", "parental")),
#   #                  structure(list(3L, c(1L, 3L), integer(0)),
#   #                            class = c("bn", "parental")),
#   #                  structure(list(2:3, 3L, integer(0)),
#   #                            class = c("bn", "parental")),
#   #                  structure(list(2L, integer(0), 1:2),
#   #                            class = c("bn", "parental")),
#   #                  structure(list(2:3, integer(0), 2L),
#   #                            class = c("bn", "parental")),
#   #                  structure(list(integer(0), 1L, 1:2),
#   #                            class = c("bn", "parental")))
#   #
#   expect_that(out, is_identical_to(expected))
# })
