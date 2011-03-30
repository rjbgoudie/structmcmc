# Part of the "structural" package, http://github.com/rbtgde/structural
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structural
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Exact Posterior methods")

test_that("Two node example", {
  x1 <- as.factor(c(1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1))
  dat <- data.frame(x1 = x1, x2 = x2)

  bnspace <- enumerateBNSpace(2, allowCyclic = TRUE)
  bnspace <- filterCyclic(bnspace)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)

  pgp <- gp(post)

  # this is exp(lsmd)/sum(exp(lsmd)), plus names
  expect_that(pgp, equals(structure(c(0.326086956521739, 0.326086956521739, 0.347826086956522), .Names = c("integer(0),1", "2,integer(0)", "integer(0),integer(0)"), class = "gp")))

  pep <- ep(post)

  # just the scores for the two graphs with single edges
  expect_that(pep, equals(structure(c(0, 0.326086956521739, 0.326086956521739, 0), .Dim = c(2L,
2L), class = c("ep", "matrix"))))
})

test_that("Second Two node example", {
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  dat <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  bnspace <- enumerateBNSpace(3, allowCyclic = TRUE)
  bnspace <- filterCyclic(bnspace)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)

  pgp <- gp(post)

  # this is exp(lsmd)/sum(exp(lsmd)), plus names
  expect_that(pgp,
    equals(
      structure(c(0.00586319286714086, 0.0439739465035564, 0.0426662664617597,
0.00586319286714086, 0.00282894696042746, 0.123127050209958,
0.123127050209958, 0.0439739465035564, 0.0426662664617595, 0.0594078861689767,
0.0426662664617595, 0.0426662664617595, 0.020586156306784, 0.00586319286714086,
0.0079210514891969, 0.123127050209958, 0.0439739465035564, 0.0426662664617595,
0.0594078861689767, 0.0426662664617597, 0.0576412376589953, 0.0079210514891969,
0.00282894696042746, 0.00274482084090454, 0.00382184844359047
), .Names = c("integer(0),1,1", "integer(0),3,1", "integer(0),c(1,3),1",
"2,integer(0),1", "integer(0),integer(0),1", "integer(0),1,2",
"2,integer(0),2", "3,integer(0),2", "2:3,integer(0),2", "integer(0),integer(0),2",
"integer(0),1,1:2", "2,integer(0),1:2", "integer(0),integer(0),1:2",
"3,1,integer(0)", "integer(0),1,integer(0)", "2,3,integer(0)",
"3,3,integer(0)", "2:3,3,integer(0)", "integer(0),3,integer(0)",
"3,c(1,3),integer(0)", "integer(0),c(1,3),integer(0)", "2,integer(0),integer(0)",
"3,integer(0),integer(0)", "2:3,integer(0),integer(0)", "integer(0),integer(0),integer(0)"
), class = "gp")
  ))

  pep <- ep(post)
  expect_that(pep[1, 3], equals(0.207114234890328)) # = sum(pgp[c(1:5, 11:13)])
  expect_that(pep[1, 1], equals(0))
})