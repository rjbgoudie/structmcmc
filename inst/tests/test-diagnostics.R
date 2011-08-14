# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Diagnostics")

test_that("cumep", {
  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1))
  dat <- data.frame(x1 = x1, x2 = x2)

  nSamples <- 5
  initial <- bn(integer(0), integer(0))

  sampler1 <- BNSampler(dat, initial, prior = priorUniform(initial))
  sampler2 <- BNSampler(dat, initial, prior = priorUniform(initial))
  sink(tempfile())
  samples1 <- draw(sampler1, n = nSamples, burnin = 0)
  
  # this is a mistake,
  # but it doesn't matter here
  # but do have to work around it, due to draw() being smarter
  samples2 <- draw(sampler2, n = nSamples, burnin = 0)
  
  sink()

  mpost1 <- bnpostmcmc(sampler1, samples1)
  mpost2 <- bnpostmcmc(sampler2, samples2)

  testmpostl <- list(mpost1, mpost2)
  class(testmpostl) <- "bnpostmcmc.list"

  expected <- list(
    # cols in order 1->1, 1->2, 2->1, 2->2
    # ie in splom's as.table order
    `Sample 1` = matrix(c(
      0, 0, 0, 0, 0,
      0/1, 1/2, 1/3, 1/4, 1/5,
      1/1, 1/2, 2/3, 2/4, 2/5,
      0, 0, 0, 0, 0
    ), ncol = 4, nrow = 5)
    ,
    # cols in order 1->1, 1->2, 2->1, 2->2
    `Sample 2` = matrix(c(
      0, 0, 0, 0, 0,
      1, 0.5, 2/3, 0.5, 3/5,
      0, 0, 0, 1/4, 1/5,
      0, 0, 0, 0, 0
    ), ncol = 4, nrow = 5)
  )
  attr(expected, "lengthOfRuns") <- 5
  attr(expected, "nbin") <- 5
  attr(expected, "numberOfNodes") <- 2
  attr(expected, "numberOfRuns") <- 2
  attr(expected, "type") <- "bnpostmcmc.list"
  class(expected) <- "epmx"
  attr(expected, "function") <- "cum"

  expect_that(
    cumep(testmpostl, method = "offline", nbin = 5),
    equals(expected)
  )

  # Moving window
  expected <- list(
    # cols in order 1->1, 1->2, 2->1, 2->2
    # ie in splom's as.table order
    `Sample 1` = matrix(c(
        0,   0,   0,
      1/3, 1/3, 0/3,
      2/3, 1/3, 1/3,
        0,   0,   0
    ), ncol = 4, nrow = 3)
    ,
    # cols in order 1->1, 1->2, 2->1, 2->2
    `Sample 2` = matrix(c(
        0,   0,   0,
      2/3, 1/3, 2/3,
      0/3, 1/3, 1/3,
        0,   0,   0
    ), ncol = 4, nrow = 3)
  )
  attr(expected, "lengthOfRuns") <- 5
  attr(expected, "nbin") <- 5
  attr(expected, "numberOfNodes") <- 2
  attr(expected, "numberOfRuns") <- 2
  attr(expected, "type") <- "bnpostmcmc.list"
  class(expected) <- "epmx"
  attr(expected, "function") <- "mw"

  if (require(zoo)){
    expect_that(
      mwep(testmpostl, method = "offline", window = 3, nbin = 5),
      equals(expected)
    )
  }

  # More
  expect_that(ep(testmpostl), equals(
    structure(list(
      structure(c(0, 0.4, 0.2, 0), .Dim = c(2L, 2L),
                class = c("ep", "matrix")),
      structure(c(0, 0.2, 0.6, 0), .Dim = c(2L, 2L),
                class = c("ep", "matrix"))),
      class = "ep.list"
    )
  ))

  expect_that(ep(testmpostl, nbin = 2, method = "flatten"), throws_error())

  expect_that(ep(testmpostl, nbin = 5, method = "flatten"), equals(
    structure(list(
      structure(list(
        structure(c(0, 1, 0, 0), .Dim = c(2L, 2L), class = c("ep", "matrix")),
        structure(c(0, 0, 1, 0), .Dim = c(2L, 2L), class = c("ep", "matrix")),
        structure(c(0, 1, 0, 0), .Dim = c(2L, 2L), class = c("ep", "matrix")),
        structure(c(0, 0, 0, 0), .Dim = c(2L, 2L), class = c("ep", "matrix")),
        structure(c(0, 0, 0, 0), .Dim = c(2L, 2L), class = c("ep", "matrix"))
      ), class = "ep.list"),
      structure(list(
        structure(c(0, 0, 1, 0), .Dim = c(2L, 2L), class = c("ep", "matrix")),
        structure(c(0, 0, 0, 0), .Dim = c(2L, 2L), class = c("ep", "matrix")),
        structure(c(0, 0, 1, 0), .Dim = c(2L, 2L), class = c("ep", "matrix")),
        structure(c(0, 1, 0, 0), .Dim = c(2L, 2L), class = c("ep", "matrix")),
        structure(c(0, 0, 1, 0), .Dim = c(2L, 2L), class = c("ep", "matrix"))
      ), class = "ep.list")
    ), class = "ep.list")
  ))

  # # Plotting tests
  #     serialised <- 'structure(list(structure(list(
  #     samples = structure(list(structure(list(
  #     2L, integer(0)), class = c("bn", "parental")), structure(list(
  #     integer(0), 1L), class = c("bn", "parental")), structure(list(
  #     2L, integer(0)), class = c("bn", "parental")), structure(list(
  #     integer(0), integer(0)), class = c("bn", "parental")), structure(list(
  #     integer(0), integer(0)), class = c("bn", "parental"))), class = c("mcmcbn",
  # "bn.list", "parental.list")), tabulated = structure(c(1L, 2L,
  # 2L), .Dim = 3L, .Dimnames = list(c("integer(0),1", "2,integer(0)",
  # "integer(0),integer(0)")), class = "table"), data = structure(list(
  #     x1 = structure(c(2L, 2L, 1L, 2L), .Label = c("0", "1"), class = "factor"),
  #     x2 = structure(c(1L, 2L, 1L, 2L), .Label = c("0", "1"),
  #     class = "factor")), .Names = c("x1",
  # "x2"), row.names = c(NA, -4L), class = "data.frame")), .Names = c("samples",
  # "tabulated", "data"), class = "bnpostmcmc"), structure(list(
  #   samples = structure(list(
  #     structure(list(integer(0), 1L), class = c("bn", "parental"
  #     )), structure(list(integer(0), integer(0)), class = c("bn",
  #     "parental")), structure(list(integer(0), 1L), class = c("bn",
  #     "parental")), structure(list(2L, integer(0)), class = c("bn",
  #     "parental")), structure(list(integer(0), 1L), class = c("bn",
  #     "parental"))), class = c("mcmcbn", "bn.list", "parental.list"
  # )), tabulated = structure(c(1L, 1L, 3L), .Dim = 3L, .Dimnames = list(
  #     c("2,integer(0)", "integer(0),integer(0)", "integer(0),1"
  #     )), class = "table"), data = structure(list(x1 = structure(c(2L,
  # 2L, 1L, 2L), .Label = c("0", "1"), class = "factor"), x2 = structure(c(1L,
  # 2L, 1L, 2L), .Label = c("0", "1"), class = "factor")), .Names = c("x1",
  # "x2"), row.names = c(NA, -4L), class = "data.frame")), .Names = c("samples",
  # "tabulated", "data"), class = "bnpostmcmc")), class = "bnpostmcmc.list")'
  #
  #   plots <- c(paste("print(splom(cumep(", serialised, ", nbin = 5)))"))
  #
  #   controlfile <- system.file("tests", "data", "diagnostics-plot-test",
  #                              package = "structmcmc")
  #   controlfn <- function(){
  #     # run this function to generate the controls
  #     # first cd into structmcmc/tests
  #     # then run this function
  #     currentwd <- getwd()
  #     #setwd("structmcmc/tests")
  #     library(graphicsQC)
  #     set.seed(301)
  #     plotcontrol <- plotExpr(plots,
  #                             path     = controlfile,
  #                             clear    = TRUE,
  # #                           filetype = c("pdf", "png"),
  #                             filetype = "png",
  #                             prefix   = "control")
  #     #setwd(currentwd)
  #   }
  #
  #   if (require(graphicsQC) & R.version$os == "darwin9.8.0"){
  #     testfile <- system.file("tests", "data", "diagnostics-plot-test",
  #                             package = "structmcmc")
  #     # generate test data
  #     set.seed(301)
  #     plottest <- plotExpr(plots,
  #                          path     = testfile,
  #                          clear    = TRUE,
  #   #                           filetype = c("pdf", "png"),
  #                               filetype = "png",
  #                          prefix   = "test")
  #
  #     # compare test data to the controls
  #     sink(tempfile())
  #     res <- compare(test = plottest, control = controlfile)
  #     sink()
  #
  #     # check that tests and controls are identical
  #
  #     for (i in seq_along(plots)){
  #       expect_that(res$results$png[[i]]$result, is_identical_to("identical"))
  #     }
  #   }
})

test_that("Fast cumep equals slow cumep", {
  set.seed(9501)
  dat <- data.frame(x1 = as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0)),
                    x2 = as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0)),
                    x3 = as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0)))

  mcmc <- posterior(data = dat, method = "gibbs", verbose = F,
                    nSamples = 1000, nBurnin = 0)

  epfast <- cumep(bnpostmcmc.list(mcmc))
  epslow <- cumep(bnpostmcmc.list(mcmc), method = "offline", nbin = 1)

  expect_equal(epfast[[1]], epslow[[1]])
})

test_that("cumtvd", {
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
  # this is a mistake,
  # but it doesn't matter here
  # but do have to work around it, due to draw() being smarter
  samples2 <- draw(sampler1, nSamples)
  
  sink()

  mpost1 <- bnpostmcmc(sampler1, samples1)
  mpost2 <- bnpostmcmc(sampler2, samples2)

  testmpostl <- list(mpost1, mpost2)
  class(testmpostl) <- "bnpostmcmc.list"

  bnspace <- enumerateBNSpace(2, allowCyclic = TRUE)
  bnspace <- filterCyclic(bnspace)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  post <- bnpost(bnspace = bnspace, data = dat, logScore = lsmd)

  # integer(0),1          2,integer(0) integer(0),integer(0) 
  #    0.3260870             0.3260870             0.3478261 
  pgp <- gp(post)

  expect_that(
    cumtvd(pgp, testmpostl, nbin = 5)$cumtvds[[2]],
    equals(
      unname(c(
      # these computed by hand
      # because we know the samples (due to set.seed)
      1 - pgp[3] + pgp[2] + pgp[1],
      0.5 - pgp[3] + pgp[2] + 0.5 - pgp[1],
      abs(2/3 - pgp[3]) + abs(0 - pgp[2]) + abs(1/3 - pgp[1]),
      abs(2/4 - pgp[3]) + abs(1/4 - pgp[2]) + abs(1/4 - pgp[1]),
      abs(3/5 - pgp[3]) + abs(1/5 - pgp[2]) + abs(1/5 - pgp[1])
      ))
    )
  )
})

test_that("cumep with sampler", {
  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1))
  dat <- data.frame(x1 = x1, x2 = x2)

  nSamples <- 5
  initial <- bn(integer(0), integer(0))

  sampler1 <- BNSampler(dat, initial, prior = priorUniform(initial))
  sampler2 <- BNSampler(dat, initial, prior = priorUniform(initial))
  sink(tempfile())
  samples1 <- draw(sampler1, n = nSamples, burnin = 0)

  # this is a mistake,
  # but it doesn't matter here
  # but do have to work around it, due to draw() being smarter
  samples2 <- draw(sampler2, n = nSamples, burnin = 0)

  sink()

  mpost1 <- bnpostmcmc(sampler1, samples1)
  mpost2 <- bnpostmcmc(sampler2, samples2)

  testmpostl <- list(mpost1, mpost2)
  class(testmpostl) <- "bnpostmcmc.list"

  expect_identical(
    cumep(testmpostl, method = "offline", nbin = 1)[[1]],
    cumep(samplers(sampler1, sampler2))[[1]]
  )
})

test_that("cumep with sampler", {
  set.seed(9501)
  dat <- data.frame(x1 = as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0)),
                    x2 = as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0)),
                    x3 = as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0)))

  mcmc1 <- posterior(data = dat, method = "mc3", verbose = F,
                     nSamples = 2000, nBurnin = 0)
  mcmc2 <- posterior(data = dat, method = "mc3", verbose = F,
                     nSamples = 2000, nBurnin = 0)

  subset <- c(2, 3)
  splom(cumep(bnpostmcmc.list(mcmc1, mcmc2)), subset = subset)
})

test_that("gelman", {
  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1))
  dat <- data.frame(x1 = x1, x2 = x2)

  nSamples <- 5000
  initial <- bn(integer(0), integer(0))

  sampler1 <- BNSampler(dat, initial, prior = priorUniform(initial))
  sampler2 <- BNSampler(dat, initial, prior = priorUniform(initial))

  samples1 <- draw(sampler1, n = nSamples, burnin = 0, verbose = F)
  samples2 <- draw(sampler2, n = nSamples, burnin = 0, verbose = F)

  expected <- c(0.99980265162473, 0.999812178389246)
  expect_equal(gelman(samplers(sampler1, sampler2)), expected)

  manual1 <- sapply(samples1, nEdges)
  manual2 <- sapply(samples2, nEdges)

  expected <- gelman.diag(mcmc.list(mcmc(manual1), mcmc(manual2)))
  expect_equal(gelman(samplers(sampler1, sampler2)), expected)

  expected <- gelman.diag(mcmc.list(mcmc(log(manual1)), mcmc(log(manual2))))
  actual <- gelman(samplers(sampler1, sampler2), transform = log)
  expect_equal(actual, expected)
})
