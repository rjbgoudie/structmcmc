# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("Statistics")

test_that("nEdges", {
  #set.seed(7101)
  set.seed(5141)
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  theData <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  fam <- enumerateBNSpace(3)
  scores <- logScoreMultDir(fam, theData)

  priors <- rep(1/25, 25)
  scores <- scores - max(scores)
  expected <- exp(scores)*priors/sum(exp(scores)*priors)

  numberOfBurnIn <- 1000
  numberOfSamples <- 60000

  expectedTable <- data.frame(expected = expected * numberOfSamples)
  row.names(expectedTable) <- lapply(fam, function(network){
    paste(network, sep = "", collapse = ",")
  })

  empty <- list(c(),c(),c())

  initial <- bn(integer(0), integer(0), integer(0))
  sampler <- BNGibbsSampler(data             = theData,
                            initial          = initial,
                            localPriors      = priorUniform(initial),
                            maxNumberParents = 2)

  samples <- draw(sampler,
                  numberOfSamples,
                  burnin = numberOfBurnIn,
                  verbose = F)
  expect_equal(statistics(sampler, "nEdges"), sapply(samples, nEdges))
})
