# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Multinomial logScoreMultDir() test")

test_that("Cooper & Herskovits 1992", {
  sink(tempfile())
  require(deal)
  sink()

  ## Test 1
  # Cooper & Herskovits 1992
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  data <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  deal <- function(){

    net <- network(data)
    sink(tempfile())
    jp <- jointprior(net)
    jp$jointalpha[1,1,1] <- 1/4
    jp$jointalpha[1,1,0] <- 1/4
    jp$jointalpha[1,1,2] <- 1/4
    jp$jointalpha[1,2,1] <- 1/4
    jp$jointalpha[1,2,2] <- 1/4
    jp$jointalpha[2,2,2] <- 1/4
    jp$jointalpha[2,1,1] <- 1/4
    jp$jointalpha[2,1,2] <- 1/4
    jp$jointalpha[2,2,1] <- 1/4
    allnetworks <- networkfamily(data,net,jp)
    allnetworks <- nwfsort(allnetworks$nw)
    sink()
    allnetworks[[1]]$score
  }

  code <- function() logScoreMultDir(bn(integer(0), 1, 2), data)

  expect_that(code(), equals(deal()))
})

test_that("Hand-computed", {
  # test against paper answer

  # this was drawn using threeNodeSampler in simple.ancst.sampl.R
  data <- structure(list(as.factor.a. = structure(c(1L, 1L, 2L, 1L, 1L,
  2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L), .Label = c("0",
  "1"), class = "factor"), as.factor.b. = structure(c(1L, 2L, 2L,
  2L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L,
  1L), .Label = c("0", "1"), class = "factor"), as.factor.c. = structure(c(1L,
  1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L,
  1L, 2L, 2L), .Label = c("0", "1"), class = "factor")), .Names = c("as.factor.a.",
  "as.factor.b.", "as.factor.c."), row.names = c(NA, -20L), class = "data.frame")

  hand <- function(){
    N_1_1_0 = 7
    N_1_1_1 = 13
    N_prime_1_1_0 = 1
    N_prime_1_1_1 = 1

    N_2_1_0 = 8
    N_2_1_1 = 12
    N_prime_2_1_0 = 1
    N_prime_2_1_1 = 1

    N_3_00_0 = 5
    N_3_00_1 = 0
    N_3_01_0 = 1
    N_3_01_1 = 7
    N_3_10_0 = 0
    N_3_10_1 = 3
    N_3_11_0 = 1
    N_3_11_1 = 3

    N_prime_3_00_0 = 1/4
    N_prime_3_00_1 = 1/4
    N_prime_3_01_0 = 1/4
    N_prime_3_01_1 = 1/4
    N_prime_3_10_0 = 1/4
    N_prime_3_10_1 = 1/4
    N_prime_3_11_0 = 1/4
    N_prime_3_11_1 = 1/4

    lastproduct1 <- (lgamma(N_prime_1_1_0 + N_1_1_0) - lgamma(N_prime_1_1_0)) + (lgamma(N_prime_1_1_1 + N_1_1_1) - lgamma(N_prime_1_1_1))

    lastproduct2 <- (lgamma(N_prime_2_1_0 + N_2_1_0) - lgamma(N_prime_2_1_0)) + (lgamma(N_prime_2_1_1 + N_2_1_1) - lgamma(N_prime_2_1_1))

    lastproduct3 <- (lgamma(N_prime_3_00_0 + N_3_00_0) - lgamma(N_prime_3_00_0)) + (lgamma(N_prime_3_01_0 + N_3_01_0) - lgamma(N_prime_3_01_0)) + (lgamma(N_prime_3_10_0 + N_3_10_0) - lgamma(N_prime_3_10_0)) + (lgamma(N_prime_3_11_0 + N_3_11_0) - lgamma(N_prime_3_11_0)) + (lgamma(N_prime_3_00_1 + N_3_00_1) - lgamma(N_prime_3_00_1)) + (lgamma(N_prime_3_01_1 + N_3_01_1) - lgamma(N_prime_3_01_1)) + (lgamma(N_prime_3_10_1 + N_3_10_1) - lgamma(N_prime_3_10_1)) + (lgamma(N_prime_3_11_1 + N_3_11_1) - lgamma(N_prime_3_11_1))

    middleproduct1 <- lgamma(N_prime_1_1_0 + N_prime_1_1_1) - lgamma(N_prime_1_1_0 + N_prime_1_1_1 + N_1_1_0 + N_1_1_1)
    middleproduct2 <- lgamma(N_prime_2_1_0 + N_prime_2_1_1) - lgamma(N_prime_2_1_0 + N_prime_2_1_1 + N_2_1_0 + N_2_1_1)

    middleproduct3 <- lgamma(N_prime_3_00_0 + N_prime_3_00_1) - lgamma(N_prime_3_00_0 + N_prime_3_00_1 + N_3_00_0 + N_3_00_1) + lgamma(N_prime_3_01_0 + N_prime_3_01_1) - lgamma(N_prime_3_01_0 + N_prime_3_01_1 + N_3_01_0 + N_3_01_1) + lgamma(N_prime_3_10_0 + N_prime_3_10_1) - lgamma(N_prime_3_10_0 + N_prime_3_10_1 + N_3_10_0 + N_3_10_1) + lgamma(N_prime_3_11_0 + N_prime_3_11_1) - lgamma(N_prime_3_11_0 + N_prime_3_11_1 + N_3_11_0 + N_3_11_1)

    loganswer <- lastproduct1 + lastproduct2 + lastproduct3 + middleproduct1 + middleproduct2 + middleproduct3

    loganswer
  }

  code <- function() logScoreMultDir(bn(integer(0), integer(0), c(1, 2)), data)

  expect_that(code(), equals(hand()))
})

test_that("Test 3", {
  # this was drawn using threeNodeSampler in simple.ancst.sampl.R
  data <- structure(list(as.factor.a. = structure(c(1L, 1L, 2L, 1L, 1L,
  2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L), .Label = c("0",
  "1"), class = "factor"), as.factor.b. = structure(c(1L, 2L, 2L,
  2L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L,
  1L), .Label = c("0", "1"), class = "factor"), as.factor.c. = structure(c(1L,
  1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L,
  1L, 2L, 2L), .Label = c("0", "1"), class = "factor")), .Names = c("as.factor.a.",
  "as.factor.b.", "as.factor.c."), row.names = c(NA, -20L), class = "data.frame")

  deal <- function(){
    sink(tempfile())
    library(deal)
    net <- network(data)
    jp <- jointprior(net)
    jp$jointalpha[1,1,1] <- 1/4
    jp$jointalpha[1,1,0] <- 1/4
    jp$jointalpha[1,1,2] <- 1/4
    jp$jointalpha[1,2,1] <- 1/4
    jp$jointalpha[1,2,2] <- 1/4
    jp$jointalpha[2,2,2] <- 1/4
    jp$jointalpha[2,1,1] <- 1/4
    jp$jointalpha[2,1,2] <- 1/4
    jp$jointalpha[2,2,1] <- 1/4
    allnetworks <- networkfamily(data,net,jp)
    allnetworks <- nwfsort(allnetworks$nw)
    sink()
    allnetworks[[1]]$score
  }

  code <- function() logScoreMultDir(bn(integer(0), integer(0), c(1, 2)), data)

  expect_that(code(), equals(code()))
})

test_that("0 node network", {
  expect_that(logScoreMultDir(bn(), data.frame()), equals(0))
  expect_that(logScoreMultDir(bn(), data.frame(c(1,2,3), c(3,2,1))), throws_error())
  expect_that(logScoreMultDir(bn(), data.frame(factor(c(1,2,3)), factor(c(3,2,1)))), equals(0))
})

test_that("1 node network", {
  testdata <- data.frame(factor(c(rep(0, 10), rep(1, 20))))
  score <- logScoreMultDir(bn(integer(0)), testdata)

  byhand_gamma <- log((gamma(2)/gamma(10 + 20 + 2)) * (gamma(10 + 1)/gamma(1)) * (gamma(20 + 1)/gamma(1)))
  byhand_lgamma <- (lgamma(2) - lgamma(10 + 20 + 2)) + (lgamma(10 + 1) - lgamma(1)) + (lgamma(20 + 1) - lgamma(1))

  expect_that(score, equals(byhand_gamma))
  expect_that(score, equals(byhand_lgamma))
})

test_that("1 node network with 0 count", {
  testdata <- data.frame(factor(c(rep(1, 30)), levels = c(0, 1)))
  score <- logScoreMultDir(bn(integer(0)), testdata)

  byhand_gamma <- log((gamma(2)/gamma(0 + 30 + 2)) * (gamma(0 + 1)/gamma(1)) * (gamma(30 + 1)/gamma(1)))
  byhand_lgamma <- (lgamma(2) - lgamma(0 + 30 + 2)) + (lgamma(0 + 1) - lgamma(1)) + (lgamma(30 + 1) - lgamma(1))

  expect_that(score, equals(byhand_gamma))
  expect_that(score, equals(byhand_lgamma))
})

test_that("1 node network with no obs", {
  testdata <- data.frame(factor(integer(0), levels = c(0, 1)))
  score <- logScoreMultDir(bn(integer(0)), testdata)

  byhand_gamma <- log((gamma(2)/gamma(0 + 0 + 2)) * (gamma(0 + 1)/gamma(1)) * (gamma(0 + 1)/gamma(1)))
  byhand_lgamma <- (lgamma(2) - lgamma(0 + 0 + 2)) + (lgamma(0 + 1) - lgamma(1)) + (lgamma(0 + 1) - lgamma(1))

  expect_that(score, equals(byhand_gamma))
  expect_that(score, equals(byhand_lgamma))
})

test_that("1 node network with no obs, 1 level", {
  testdata <- data.frame(factor(integer(0), levels = c(0)))
  score <- logScoreMultDir(bn(integer(0)), testdata)

  byhand_gamma <- log((gamma(1)/gamma(0 + 0 + 1)) * (gamma(0 + 1)/gamma(1)) * (gamma(0 + 1)/gamma(1)))
  byhand_lgamma <- (lgamma(1) - lgamma(0 + 0 + 1)) + (lgamma(0 + 1) - lgamma(1)) + (lgamma(0 + 1) - lgamma(1))

  expect_that(score, equals(byhand_gamma))
  expect_that(score, equals(byhand_lgamma))
})

test_that("1 node network with no obs, 1 level", {
  testdata <- data.frame(factor(c(1:10)))
  score <- logScoreMultDir(bn(integer(0)), testdata)

  firstfactor_gamma <- function(Nij, Nij_prime) gamma(Nij_prime)/gamma(Nij + Nij_prime)
  firstfactor_lgamma <- function(Nij, Nij_prime) lgamma(Nij_prime) - lgamma(Nij + Nij_prime)

  lastfactor_gamma <- function(Nijk, Nijk_prime) gamma(Nijk + Nijk_prime)/gamma(Nijk_prime)
  lastfactor_lgamma <- function(Nijk, Nijk_prime) lgamma(Nijk + Nijk_prime) - lgamma(Nijk_prime)

  byhand_gamma <- log(firstfactor_gamma(10, 10) * prod(rep(lastfactor_gamma(1, 1), 10)))
  byhand_lgamma <- firstfactor_lgamma(10, 10) + sum(rep(lastfactor_lgamma(1, 1), 10))

  expect_that(score, equals(byhand_gamma))
  expect_that(score, equals(byhand_lgamma))
})

test_that("Sewell and Shah data via Heckerman A Tutorial on Learning With Bayesian Networks", {
  sufficient_statistics <- c(4, 349, 13, 64, 9, 207, 33, 72, 12, 126, 38, 54, 10, 67, 49, 43, 2, 232, 27, 84, 7, 201, 64, 95, 12, 115, 93, 92, 17, 79, 119, 59, 8, 166, 47, 91, 6, 120, 74, 110, 17, 92, 148, 100, 6, 42, 198, 73, 4, 48, 39, 57, 5, 47, 123, 90, 9, 41, 224, 65, 8, 17, 414, 54, 5, 454, 9, 44, 5, 312, 14, 47, 8, 216, 20, 35, 13, 96, 28, 24, 11, 285, 29, 61, 19, 236, 47, 88, 12, 164, 62, 85, 15, 113, 72, 50, 7, 163, 36, 72, 13, 193, 75, 90, 12, 174, 91, 100, 20, 81, 142, 77, 6, 50, 36, 58, 5, 70, 110, 76, 12, 48, 230, 81, 13, 49, 360, 98)

  levs <- expand.grid(CP = c("Yes", "No"), PE = c("Low", "High"), IQ = c("Low", "Lower Middle", "Upper Middle", "High"), SES = c("Low", "Lower Middle", "Upper Middle", "High"), gender = c("Male", "Female"))

  data <- data.frame()

  for (i in 1:nrow(levs)){
    data <- rbind(data, do.call("rbind", lapply(seq_len(sufficient_statistics[i]), function(x) levs[i, ])))
  }

  # need to implement the equivalent sample size hyperparameters
  # here, Heckerman uses equivalent sample size = 5
})

test_that("1 node network, hyperparameters = 0.9", {
  testdata <- data.frame(factor(c(rep(0, 10), rep(1, 20))))
  score <- logScoreMultDir(bn(integer(0)),
                           testdata,
                           hyperparameters = "point9")

  byhand_gamma <- log((gamma(0.9 * 2)/gamma(10 + 20 + 0.9 * 2)) * (gamma(10 + 0.9)/gamma(0.9)) * (gamma(20 + 0.9)/gamma(0.9)))
  byhand_lgamma <- (lgamma(0.9 * 2) - lgamma(10 + 20 + 0.9 * 2)) + (lgamma(10 + 0.9) - lgamma(0.9)) + (lgamma(20 + 0.9) - lgamma(0.9))

  expect_that(score, equals(byhand_gamma))
  expect_that(score, equals(byhand_lgamma))
})

test_that("Hand-computed, hyperparameters = 0.9", {
  # test against paper answer

  # this was drawn using threeNodeSampler in simple.ancst.sampl.R
  data <- structure(list(as.factor.a. = structure(c(1L, 1L, 2L, 1L, 1L,
  2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L), .Label = c("0",
  "1"), class = "factor"), as.factor.b. = structure(c(1L, 2L, 2L,
  2L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L,
  1L), .Label = c("0", "1"), class = "factor"), as.factor.c. = structure(c(1L,
  1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L,
  1L, 2L, 2L), .Label = c("0", "1"), class = "factor")), .Names = c("as.factor.a.",
  "as.factor.b.", "as.factor.c."), row.names = c(NA, -20L), class = "data.frame")

  hand <- function(){
    N_1_1_0 = 7
    N_1_1_1 = 13
    N_prime_1_1_0 = 0.9
    N_prime_1_1_1 = 0.9

    N_2_1_0 = 8
    N_2_1_1 = 12
    N_prime_2_1_0 = 0.9
    N_prime_2_1_1 = 0.9

    N_3_00_0 = 5
    N_3_00_1 = 0
    N_3_01_0 = 1
    N_3_01_1 = 7
    N_3_10_0 = 0
    N_3_10_1 = 3
    N_3_11_0 = 1
    N_3_11_1 = 3

    N_prime_3_00_0 = .9
    N_prime_3_00_1 = .9
    N_prime_3_01_0 = .9
    N_prime_3_01_1 = .9
    N_prime_3_10_0 = .9
    N_prime_3_10_1 = .9
    N_prime_3_11_0 = .9
    N_prime_3_11_1 = .9

    lastproduct1 <- (lgamma(N_prime_1_1_0 + N_1_1_0) - lgamma(N_prime_1_1_0)) + (lgamma(N_prime_1_1_1 + N_1_1_1) - lgamma(N_prime_1_1_1))

    lastproduct2 <- (lgamma(N_prime_2_1_0 + N_2_1_0) - lgamma(N_prime_2_1_0)) + (lgamma(N_prime_2_1_1 + N_2_1_1) - lgamma(N_prime_2_1_1))

    lastproduct3 <- (lgamma(N_prime_3_00_0 + N_3_00_0) - lgamma(N_prime_3_00_0)) + (lgamma(N_prime_3_01_0 + N_3_01_0) - lgamma(N_prime_3_01_0)) + (lgamma(N_prime_3_10_0 + N_3_10_0) - lgamma(N_prime_3_10_0)) + (lgamma(N_prime_3_11_0 + N_3_11_0) - lgamma(N_prime_3_11_0)) + (lgamma(N_prime_3_00_1 + N_3_00_1) - lgamma(N_prime_3_00_1)) + (lgamma(N_prime_3_01_1 + N_3_01_1) - lgamma(N_prime_3_01_1)) + (lgamma(N_prime_3_10_1 + N_3_10_1) - lgamma(N_prime_3_10_1)) + (lgamma(N_prime_3_11_1 + N_3_11_1) - lgamma(N_prime_3_11_1))

    middleproduct1 <- lgamma(N_prime_1_1_0 + N_prime_1_1_1) - lgamma(N_prime_1_1_0 + N_prime_1_1_1 + N_1_1_0 + N_1_1_1)
    middleproduct2 <- lgamma(N_prime_2_1_0 + N_prime_2_1_1) - lgamma(N_prime_2_1_0 + N_prime_2_1_1 + N_2_1_0 + N_2_1_1)

    middleproduct3 <- lgamma(N_prime_3_00_0 + N_prime_3_00_1) - lgamma(N_prime_3_00_0 + N_prime_3_00_1 + N_3_00_0 + N_3_00_1) + lgamma(N_prime_3_01_0 + N_prime_3_01_1) - lgamma(N_prime_3_01_0 + N_prime_3_01_1 + N_3_01_0 + N_3_01_1) + lgamma(N_prime_3_10_0 + N_prime_3_10_1) - lgamma(N_prime_3_10_0 + N_prime_3_10_1 + N_3_10_0 + N_3_10_1) + lgamma(N_prime_3_11_0 + N_prime_3_11_1) - lgamma(N_prime_3_11_0 + N_prime_3_11_1 + N_3_11_0 + N_3_11_1)

    loganswer <- lastproduct1 + lastproduct2 + lastproduct3 + middleproduct1 + middleproduct2 + middleproduct3

    loganswer
  }

  code <- function() logScoreMultDir(bn(integer(0), integer(0), c(1, 2)),
                                     data,
                                     hyperparameters = "point9")

  expect_that(code(), equals(hand()))
})
