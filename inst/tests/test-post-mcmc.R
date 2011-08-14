# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("MCMC Posterior methods")

test_that("Two node example", {
  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1))
  dat <- data.frame(x1 = x1, x2 = x2)

  nSamples <- 5
  initial <- bn(integer(0), integer(0))

  sampler <- BNSampler(dat, initial, prior = priorUniform(initial))
  sink(tempfile())
  samples <- draw(sampler, nSamples)
  sink()

  mpost <- bnpostmcmc(sampler, samples)

  mgp <- gp(mpost)

  expect_that(mgp, equals(
    structure(c(0.2, 0.4, 0.4), .Names = c("integer(0),1", "2,integer(0)",
"integer(0),integer(0)"), class = "gp")
  ))

  mep <- ep(mpost)

  expect_that(mep, equals(
    structure(c(0, 0.4, 0.2, 0), .Dim = c(2L, 2L), class = c("ep",
"matrix"))
  ))
  expect_identical(mep, ep(sampler))
})

test_that("Tabulate samples", {
  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1))
  dat <- data.frame(x1 = x1, x2 = x2)

  nSamples <- 5
  initial <- bn(integer(0), integer(0))

  sampler1 <- BNSampler(dat, initial, prior = priorUniform(initial))
  sampler2 <- BNSampler(dat, initial, prior = priorUniform(initial))
  sink(tempfile())
  samples1 <- draw(sampler1, nSamples)
  samples2 <- draw(sampler1, nSamples) # this is a mistake, but it doesn't matter here
  sink()

  expect_that(pltabulate(samples1, sort = T), equals(
    structure(c(1L, 2L, 2L), .Dim = 3L, .Dimnames = structure(list(
    c("integer(0),1", "2,integer(0)", "integer(0),integer(0)"
    ))), class = "table")
  ))

  expect_that(pltabulate(samples1, pretty = T, sort = T), equals(
    structure(c(1L, 2L, 2L), .Dim = 3L, .Dimnames = list(c("[][1]",
"[][]", "[2][]")), class = "table")
  ))
})

test_that("map and top", {
  set.seed(1331)
  x1 <- as.factor(c(1, 1, 0, 1, 2, 1, 2, 3))
  x2 <- as.factor(c(0, 1, 0, 1, 1, 1, 1, 2))
  x3 <- as.factor(c(2, 0, 1, 1, 0, 0, 1, 2))
  dat <- data.frame(x1 = x1, x2 = x2, x3 = x3)

  nSamples <- 20
  initial <- bn(integer(0), integer(0), integer(0))

  sampler <- BNSampler(dat, initial, prior = priorUniform(initial))
  sink(tempfile())
  samples <- draw(sampler, nSamples)
  sink()

  mpost <- bnpostmcmc(sampler = sampler, samples = samples)

  expect_that(map(mpost), is_identical_to(bn(2L, 3L, integer(0))))

  expect_that(top(mpost, head = 2),
    is_identical_to(bn.list(bn(2L, 3L, integer(0)),
                            bn(2L, integer(0), integer(0)))))

  # these expectations differ between platforms due to differences in
  # collation order.
  # http://markmail.org/message/4q23tyzh3hzp744w
  # sort yields different results on OS X (PR#14163) - Prof Brian Ripley
  # if (R.version$os == "darwin9.8.0"){
    expect_that(logScoreMultDir(mpost, sampler, dat),
                equals(c(-29.4539925134359,
                         -30.0375534255437,
                         -31.6857311662937,
                         -29.4539925134359,
                         -31.9150536235599,
                         -29.3665002516089,
                         -29.9500611637167,
                         -29.8073335035549,
                         -30.9045731084110)))
  # } else {
  #   expect_that(logScoreMultDir(mpost, sampler, dat),
  #               equals(c(-29.4539925134359,
  #                        -30.0375534255437,
  #                        -31.9150536235599,
  #                        -31.6857311662937,
  #                        -29.4539925134359,
  #                        -29.3665002516089,
  #                        -29.9500611637167,
  #                        -29.8073335035549,
  #                        -30.9045731084110)))
  # }

  expect_that(logScoreMultDir(mpost, sampler, dat,
                              sort.by = "logScoreMultDir"),
              equals(c(-29.3665002516089,
                       -29.4539925134359,
                       -29.4539925134359,
                       -29.8073335035549,
                       -29.9500611637167,
                       -30.0375534255437,
                       -30.9045731084110,
                       -31.6857311662937,
                       -31.9150536235599)))

  # http://markmail.org/message/4q23tyzh3hzp744w
  # sort yields different results on OS X (PR#14163) - Prof Brian Ripley
  # if (R.version$os == "darwin9.8.0"){
    expect_that(logScoreMultDir(mpost, sampler, dat,
                              use.names = T),
                equals(structure(c(-29.4539925134359,
                                   -30.0375534255437,
                                   -31.6857311662937,
                                   -29.4539925134359,
                                   -31.9150536235599,
                                   -29.3665002516089,
                                   -29.9500611637167,
                                   -29.8073335035549,
                                   -30.9045731084110),
                        .Names = c("2,3,integer(0)",
                                   "2,integer(0),integer(0)",
                                   "2,integer(0),1",
                                   "2,integer(0),2",
                                   "2:3,integer(0),integer(0)",
                                   "integer(0),1,2",
                                   "integer(0),1,integer(0)",
                                   "integer(0),c(1,3),integer(0)",
                                   "integer(0),integer(0),integer(0)"))))
  # } else {
  #     expect_that(logScoreMultDir(mpost, sampler, dat,
  #                               use.names = T),
  #               equals(structure(c(-29.4539925134359,
  #                                  -30.0375534255437,
  #                                  -31.9150536235599,
  #                                  -31.6857311662937,
  #                                  -29.4539925134359,
  #                                  -29.3665002516089,
  #                                  -29.9500611637167,
  #                                  -29.8073335035549,
  #                                  -30.9045731084110),
  #                         .Names = c("2,3,integer(0)",
  #                                    "2,integer(0),integer(0)",
  #                                    "2:3,integer(0),integer(0)",
  #                                    "2,integer(0),1",
  #                                    "2,integer(0),2",
  #                                    "integer(0),1,2",
  #                                    "integer(0),1,integer(0)",
  #                                    "integer(0),c(1,3),integer(0)",
  #                                    "integer(0),integer(0),integer(0)"))))
  #   }

  expect_that(logScoreMultDir(mpost, sampler, dat,
                              sort.by = "logScoreMultDir", use.names = T),
              equals(structure(c(-29.3665002516089,
                                 -29.4539925134359,
                                 -29.4539925134359,
                                 -29.8073335035549,
                                 -29.9500611637167,
                                 -30.0375534255437,
                                 -30.9045731084110,
                                 -31.6857311662937,
                                 -31.9150536235599),
                               .Names = c("integer(0),1,2",
                                          "2,3,integer(0)",
                                          "2,integer(0),2",
                                          "integer(0),c(1,3),integer(0)",
                                          "integer(0),1,integer(0)",
                                          "2,integer(0),integer(0)",
                                          "integer(0),integer(0),integer(0)",
                                          "2,integer(0),1",
                                          "2:3,integer(0),integer(0)"))))

  expect_that(logScoreMultDir(mpost, sampler, dat,
                              sort.by = "logScoreMultDir",
                              head = 4,
                              use.names = T),
              equals(structure(c(-29.3665002516089,
                                 -29.4539925134359,
                                 -29.4539925134359,
                                 -29.8073335035549),
                               .Names = c("integer(0),1,2",
                                          "2,3,integer(0)",
                                          "2,integer(0),2",
                                          "integer(0),c(1,3),integer(0)"))))

  expect_that(logScoreMultDir(mpost, sampler, dat,
                              sort.by = "logScoreMultDir",
                              head = 3,
                              use.names = T),
              equals(structure(c(-29.3665002516089,
                                 -29.4539925134359,
                                 -29.4539925134359),
                               .Names = c("integer(0),1,2",
                                          "2,3,integer(0)",
                                          "2,integer(0),2"))))

  expect_that(logScoreMultDir(mpost, sampler, dat,
                              sort.by = "logScoreMultDir",
                              head = 2,
                              use.names = T),
              equals(structure(c(-29.3665002516089,
                                 -29.4539925134359,
                                 -29.4539925134359),
                               .Names = c("integer(0),1,2",
                                          "2,3,integer(0)",
                                          "2,integer(0),2"))))

  expect_that(logScoreMultDir(mpost, sampler, dat,
                              sort.by = "logScoreMultDir",
                              head = 1,
                              use.names = T),
              equals(structure(c(-29.3665002516089),
                               .Names = c("integer(0),1,2"))))

  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1))
  dat <- data.frame(x1 = x1, x2 = x2)

  nSamples <- 5
  initial <- bn(integer(0), integer(0))

  sampler <- BNSampler(dat, initial, prior = priorUniform(initial))
  sink(tempfile())
  samples <- draw(sampler, nSamples)
  sink()

  mpost <- bnpostmcmc(sampler, samples)

  expect_that(map(mpost), is_identical_to(
    bn.list(bn(2L, integer(0)),
            bn(integer(0), integer(0)))))

  expect_that(top(mpost, head = 3), is_identical_to(
    bn.list(bn(2L, integer(0)),
            bn(integer(0), integer(0)),
            bn(integer(0), 1))))
})
