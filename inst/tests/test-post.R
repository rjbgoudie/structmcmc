# Part of the "structural" package, http://github.com/rbtgde/structural
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structural
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Posterior methods")

test_that("ROC Curve plotting", {
  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1))
  dat <- data.frame(x1 = x1, x2 = x2)

  nSamples <- 5
  prior <- function(net) 1
  initial <- bn(integer(0), integer(0))

  sampler1 <- BNSampler(dat, initial, prior)
  sampler2 <- BNSampler(dat, initial, prior)
  sink(tempfile())
  samples1 <- draw(sampler1, nSamples)
  samples2 <- draw(sampler1, nSamples) # this is a mistake
                                       # but it doesn't matter here
  sink()

  mpost1 <- bnpostmcmc(sampler1, samples1)
  mpost2 <- bnpostmcmc(sampler2, samples2)

  testmpostl <- list(mpost1, mpost2)
  class(testmpostl) <- "bnpostmcmc.list"

  # IMPORTANT:
  # This is a weird hack to get around scoping issues with graphicsQC
  # Instead dput the testmpostl to serialise this.
  serialised <- 'structure(list(structure(list(
    samples = structure(list(structure(list(
    2L, integer(0)), class = c("bn", "parental")), structure(list(
    integer(0), 1L), class = c("bn", "parental")), structure(list(
    2L, integer(0)), class = c("bn", "parental")), structure(list(
    integer(0), integer(0)), class = c("bn", "parental")), structure(list(
    integer(0), integer(0)), class = c("bn", "parental"))), class = c("mcmcbn",
"bn.list", "parental.list")), tabulated = structure(c(1L, 2L,
2L), .Dim = 3L, .Dimnames = list(c("integer(0),1", "2,integer(0)",
"integer(0),integer(0)")), class = "table"), data = structure(list(
    x1 = structure(c(2L, 2L, 1L, 2L), .Label = c("0", "1"), class = "factor"),
    x2 = structure(c(1L, 2L, 1L, 2L), .Label = c("0", "1"),
    class = "factor")), .Names = c("x1",
"x2"), row.names = c(NA, -4L), class = "data.frame")), .Names = c("samples",
"tabulated", "data"), class = "bnpostmcmc"), structure(list(
  samples = structure(list(
    structure(list(integer(0), 1L), class = c("bn", "parental"
    )), structure(list(integer(0), integer(0)), class = c("bn",
    "parental")), structure(list(integer(0), 1L), class = c("bn",
    "parental")), structure(list(2L, integer(0)), class = c("bn",
    "parental")), structure(list(integer(0), 1L), class = c("bn",
    "parental"))), class = c("mcmcbn", "bn.list", "parental.list"
)), tabulated = structure(c(1L, 1L, 3L), .Dim = 3L, .Dimnames = list(
    c("2,integer(0)", "integer(0),integer(0)", "integer(0),1"
    )), class = "table"), data = structure(list(x1 = structure(c(2L,
2L, 1L, 2L), .Label = c("0", "1"), class = "factor"), x2 = structure(c(1L,
2L, 1L, 2L), .Label = c("0", "1"), class = "factor")), .Names = c("x1",
"x2"), row.names = c(NA, -4L), class = "data.frame")), .Names = c("samples",
"tabulated", "data"), class = "bnpostmcmc")), class = "bnpostmcmc.list")'

  set.seed(231)
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0))
  x3 <- as.factor(c(1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0))
  dat <- data.frame(x1 = x1, x2 = x2, x3 = x3)

  nSamples <- 10
  prior <- function(net) 1
  initial1 <- bn(integer(0), integer(0), integer(0))
  initial2 <- bn(c(2L, 3L), 3L, integer(0))

  samplerb1 <- BNSampler(dat, initial1, prior)
  samplerb2 <- BNSampler(dat, initial2, prior)
  sink(tempfile())
  samplesb1 <- draw(samplerb1, nSamples)
  samplesb2 <- draw(samplerb2, nSamples) # this is a mistake
                                         # but it doesn't matter here
  sink()

  mpostb1 <- bnpostmcmc(samplerb1, samplesb1)
  mpostb2 <- bnpostmcmc(samplerb2, samplesb2)

  testmpostlb <- list(mpostb1, mpostb2)
  class(testmpostlb) <- "bnpostmcmc.list"

  serialisedb <- 'structure(list(structure(list(samples = structure(
    list(structure(list(
    2L, integer(0), integer(0)), class = c("bn", "parental")),
    structure(list(2:3, integer(0), integer(0)), class = c("bn",
    "parental")), structure(list(2:3, integer(0), 2L), class = c("bn",
    "parental")), structure(list(2L, integer(0), 1:2), class = c("bn",
    "parental")), structure(list(2L, integer(0), 1L), class = c("bn",
    "parental")), structure(list(2L, integer(0), 1L), class = c("bn",
    "parental")), structure(list(2L, integer(0), integer(0)), class = c("bn",
    "parental")), structure(list(integer(0), integer(0), integer(0)),
    class = c("bn",
    "parental")), structure(list(integer(0), integer(0), integer(0)),
    class = c("bn",
    "parental")), structure(list(integer(0), integer(0), integer(0)),
    class = c("bn",
    "parental"))), class = c("mcmcbn", "bn.list", "parental.list"
)), tabulated = structure(c(1L, 1L, 1L, 2L, 2L, 3L), .Dim = 6L,
.Dimnames = list(
    c("2,integer(0),1:2", "2:3,integer(0),2", "2:3,integer(0),integer(0)",
    "2,integer(0),1", "2,integer(0),integer(0)",
    "integer(0),integer(0),integer(0)"
    )), class = "table"), data = structure(list(x1 = structure(c(2L,
2L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 2L), .Label = c("0",
"1"), class = "factor"), x2 = structure(c(1L, 2L, 1L, 2L, 2L,
2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L), .Label = c("0", "1"), class = "factor"),
    x3 = structure(c(2L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L,
    2L, 1L, 1L), .Label = c("0", "1"), class = "factor")), .Names = c("x1",
"x2", "x3"), row.names = c(NA, -13L), class = "data.frame")),
.Names = c("samples",
"tabulated", "data"), class = "bnpostmcmc"),
structure(list(samples = structure(list(
    structure(list(2:3, integer(0), 2L), class = c("bn", "parental"
    )), structure(list(2L, integer(0), 2L), class = c("bn", "parental"
    )), structure(list(2L, integer(0), integer(0)), class = c("bn",
    "parental")), structure(list(2L, integer(0), 1L), class = c("bn",
    "parental")), structure(list(integer(0), integer(0), 1L), class = c("bn",
    "parental")), structure(list(integer(0), integer(0), integer(0)),
    class = c("bn",
    "parental")), structure(list(integer(0), integer(0), integer(0)),
    class = c("bn",
    "parental")), structure(list(integer(0), integer(0), 1L), class = c("bn",
    "parental")), structure(list(integer(0), integer(0), 1L), class = c("bn",
    "parental")), structure(list(integer(0), integer(0), 1L), class = c("bn",
    "parental"))), class = c("mcmcbn", "bn.list", "parental.list"
)), tabulated = structure(c(1L, 1L, 1L, 1L, 2L, 4L), .Dim = 6L,
.Dimnames = list(
    c("2,integer(0),1", "2,integer(0),2", "2,integer(0),integer(0)",
    "2:3,integer(0),2", "integer(0),integer(0),integer(0)",
    "integer(0),integer(0),1"
    )), class = "table"), data = structure(list(x1 = structure(c(2L,
2L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 2L), .Label = c("0",
"1"), class = "factor"), x2 = structure(c(1L, 2L, 1L, 2L, 2L,
2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L), .Label = c("0", "1"), class = "factor"),
    x3 = structure(c(2L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L,
    2L, 1L, 1L), .Label = c("0", "1"), class = "factor")), .Names = c("x1",
"x2", "x3"), row.names = c(NA, -13L), class = "data.frame")),
.Names = c("samples",
"tabulated", "data"), class = "bnpostmcmc")), class = "bnpostmcmc.list")'

  plots <- c(paste("print(rocplot(true = bn(integer(0), 1L), mcmcs = ",
                 serialised,
                 "))"),
             paste("print(rocplot(true = bn(2L, integer(0), integer(0)),",
                   " mcmcs = ",
                   serialisedb,
                   ", map = bn.list(bn(integer(0), 1L, 1L))))"))

  # for example 2, the 'MAP estimate' is of course not accurate,
  # it is just an example for testing.

  controlfile <- file.path("", "Volumes", "Buster", "library",
                           "structural", "inst", "test-data",
                           "struct-dag-inf2-plot-control")

  controlfn <- function(){
    # run this function to generate the controls
    # first cd into struct-dag-inf2/tests
    # then run this function
    currentwd <- getwd()
    setwd("/Volumes/buster.stats.warwick.ac.uk/library/struct-dag-inf2/tests")
    library(graphicsQC)
    set.seed(301)
    plotcontrol <- plotExpr(plots,
                            path     = controlfile,
                            clear    = TRUE,
#                           filetype = c("pdf", "png"),
                            filetype = "png",
                            prefix   = "control")
    setwd(currentwd)
  }

  if (require(graphicsQC) && R.version$os == "darwin9.8.0"){
    testfile <- file.path("", "Volumes", "Buster", "library",
                             "structural", "inst", "test-data",
                             "struct-dag-inf2-plot-test")
    # generate test data
    set.seed(301)
    plottest <- plotExpr(plots,
                         path     = testfile,
                         clear    = TRUE,
  #                           filetype = c("pdf", "png"),
                              filetype = "png",
                         prefix   = "test")

    # compare test data to the controls
    sink(tempfile())
    res <- compare(test = plottest, control = controlfile)
    sink()

    # check that tests and controls are identical

    for (i in seq_along(plots)){
      expect_that(res$results$png[[i]]$result, is_identical_to("identical"))
    }

    expect_that(print(rocplot(true = bn(integer(0), 1),
                      map = bn(integer(0), integer(0)))),
                throws_error())
  }
})

test_that("ep.parental.list", {
  testpl <- parental.list(parental(2L, integer(0)),
                          parental(integer(0), 1L),
                          parental(integer(0), integer(0)),
                          parental(2L, integer(0)),
                          parental(2L, integer(0)),
                          parental(2L, integer(0)))

  expected <- matrix(c(0, 1/6,
                       4/6, 0),
                     nrow  = 2,
                     ncol  = 2,
                     byrow = T)
  class(expected) <- c("ep", "matrix")

  expect_that(ep(testpl, nbin = 1), is_identical_to(expected))

  # using start and end
  expected <- matrix(c(0, 1/2,
                       0, 0),
                     nrow  = 2,
                     ncol  = 2,
                     byrow = T)
  class(expected) <- c("ep", "matrix")
  expect_that(ep(testpl, nbin = 1, start = 2, end = 3),
              is_identical_to(expected))
  # errors
  expect_that(ep(testpl, nbin = 1, start = 2, end = 10),
              throws_error())
  expect_that(ep(testpl, nbin = 1, start = 0, end = 2),
              throws_error())
  expect_that(ep(testpl, nbin = 1.3),
              throws_error())
  expect_that(ep(testpl, start = 1.3),
              throws_error())
  expect_that(ep(testpl, end = 1.3),
              throws_error())

  expected1 <- matrix(c(0, 1/2,
                       1/2, 0),
                     nrow  = 2,
                     ncol  = 2,
                     byrow = T)
  class(expected1) <- c("ep", "matrix")
  expected2 <- matrix(c(0, 0,
                       1/2, 0),
                     nrow  = 2,
                     ncol  = 2,
                     byrow = T)
  class(expected2) <- c("ep", "matrix")
  expected3 <- matrix(c(0, 0,
                       1, 0),
                     nrow  = 2,
                     ncol  = 2,
                     byrow = T)
  class(expected3) <- c("ep", "matrix")

  expected <- ep.list(expected1, expected2, expected3)

  expect_that(ep(testpl, nbin = 3), is_identical_to(expected))
})

test_that("ep.table", {
  tab <- c(`integer(0),c(1L,3L),integer(0)` = 3,
           `integer(0),c(1L),integer(0)` = 6,
           `integer(0),integer(0),integer(0)` = 11)
  tab <- as.table(tab)

  sink(tempfile())
  expected <- matrix(c(0, 9/20, 0,
                       0,    0, 0,
                       0, 3/20, 0), nrow = 3, ncol = 3, byrow = TRUE)
  class(expected) <- c("ep", "matrix")
  expect_that(ep(tab),
              is_identical_to(expected))

  addEdge2to1 <- function(pl){
    stopifnot("parental.list" %in% class(pl))
    res <- lapply(pl, function(p){
      p[[1]] <- 2L
      p
    })
    class(res) <- "parental.list"
    res
  }

  expected <- matrix(c(0, 9/20, 0,
                       1,    0, 0,
                       0, 3/20, 0), nrow = 3, ncol = 3, byrow = TRUE)
  class(expected) <- c("ep", "matrix")
  expect_that(ep(tab, FUN = addEdge2to1),
              is_identical_to(expected))

  # error
  addEdge2to1NotReturningParentalList <- function(pl){
    stopifnot("parental.list" %in% class(pl))
    lapply(pl, function(p){
      p[[1]] <- 2L
      p
    })
  }

  expected <- matrix(c(0, 9/20, 0,
                       1,    0, 0,
                       0, 3/20, 0), nrow = 3, ncol = 3, byrow = TRUE)
  class(expected) <- c("ep", "matrix")
  expect_that(ep(tab, FUN = addEdge2to1NotReturningParentalList),
              throws_error())
  sink()
})

test_that("parentalFromEPThreshold", {
  testdata <- matrix(c(0, 1/6,
                       4/6, 0),
                     nrow  = 2,
                     ncol  = 2,
                     byrow = T)
  class(testdata) <- c("ep", "matrix")

  expected1 <- parental(integer(0), integer(0))
  expected2 <- parental(2L, integer(0))
  expected3 <- parental(2L, 1L)

  expect_that(parentalFromEPThreshold(testdata, 1),
              is_identical_to(expected1))
  expect_that(parentalFromEPThreshold(testdata, 0.9),
              is_identical_to(expected1))
  expect_that(parentalFromEPThreshold(testdata, 3/6),
              is_identical_to(expected2))
  expect_that(parentalFromEPThreshold(testdata, 1/6),
              is_identical_to(expected3))
  expect_that(parentalFromEPThreshold(testdata, 0),
              is_identical_to(expected3))
})

test_that("roc with bnpost", {
  x1 <- as.factor(c(1, 1, 0, 1, 1, 1))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1))
  x3 <- as.factor(c(0, 1, 0, 0, 0, 0))
  dat <- data.frame(x1 = x1, x2 = x2, x3 = x3)

  bnspace <- enumerateBNSpace(3, allowCyclic = TRUE)
  bnspace <- filterCyclic(bnspace)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  post <- bnpost(bnspace, logScore = lsmd, data = dat)
  eppost <- ep(post)

  rocplot(true = bn(integer(0), c(1,3), integer(0)),
          eps = eppost)
  rocplot(true = bn(integer(0), integer(0), 2),
          eps = eppost)
  anep <- matrix(rnorm(9), 3, 3)
  class(anep) <- c("ep", "matrix")
  rocplot(true = bn(integer(0), integer(0), 2),
          eps = ep.list(eppost, anep))
})