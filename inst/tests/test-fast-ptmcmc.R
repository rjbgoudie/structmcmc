# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("MCMC BN PT Sampling (Fast Tests)")

test_that("2-node Bayesian Network", {
  set.seed(1424)
  x1 <- as.factor(c(1, 1, 0, 1))
  x2 <- as.factor(c(0, 1, 0, 1))
  theData <- data.frame(x1 = x1, x2 = x2)

  fam <- enumerateBNSpace(2)
  scores <- logScoreMultDir(fam, theData)

  priors <- 1
  scores <- scores - max(scores)
  expected <- exp(scores)*priors/sum(exp(scores)*priors)

  numberOfBurnIn <- 100
  numberOfSamples <- 1000

  expectedTable <- data.frame(expected = expected * numberOfSamples)
  row.names(expectedTable) <- lapply(fam, function(network){
    paste(network, sep = "", collapse = ",")})

  priorFlat <- function(network) {
    1/length(fam)
  }

  cache <- new.env(hash = T, size = 10000L)
  sampler1 <- BNSamplerPT(theData, bn(integer(0), integer(0)), priorFlat, temp = 1, cache = cache)
  sampler2 <- BNSamplerPT(theData, bn(integer(0), integer(0)), priorFlat, temp = 0.5, cache = cache)
  sampler3 <- BNSamplerPT(theData, bn(integer(0), integer(0)), priorFlat, temp = 0.1, cache = cache)

  samplers <- list(sampler1, sampler2, sampler3)
  samples <- drawPT(samplers = samplers,
                     data = theData,
                     n = 1000,
                     temperatures = c(1, 0.5, 0.1),
                     pswap = 0.2,
                     burnin = 0,
                     verbose = F)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")}))))

  expect_that(as.vector(outTable["integer(0),1"]),
              is_within(326, 35))
  expect_that(as.vector(outTable["2,integer(0)"]),
              is_within(326, 35))
  expect_that(as.vector(outTable["integer(0),integer(0)"]),
              is_within(347, 20))
})
