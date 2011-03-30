# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Cross-validation tests")

test_that("3-node Bayesian Network", {

  if (require(reshape) && require(nnet)){
    x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
    x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
    x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
    theData <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

    modal_graph <- bn(integer(0), c(1,3), integer(0))

    set.seed(3012)
    cv <- dagcv(modal_graph, 2, theData, fold = list(c(seq_len(5)), c(seq_len(5) + 5)), predictiveReplications = 1)

    set.seed(3012)
    #pt2
    sample.int(2, replace = T, prob = c(9/10, 1/10), size = 1)
    sample.int(2, replace = T, prob = c(1/6, 5/6), size = 1)
    sample.int(2, replace = T, prob = c(1/2, 1/2), size = 1)
    sample.int(2, replace = T, prob = c(1/10, 9/10), size = 2)

    sample.int(2, replace = T, size = 5)

    #pt1
    sample.int(2, replace = T, prob = c(5/6, 1/6), size = 2)
    sample.int(2, replace = T, prob = c(5/6, 1/6), size = 1)
    sample.int(2, replace = T, prob = c(1/10, 9/10), size = 2)

    expect_that(cv$overall["bn"], equals(structure(3, .Names = "bn")))
    expect_that(cv$numberIncorrectByFold, equals(c(2,1)))
  }
})