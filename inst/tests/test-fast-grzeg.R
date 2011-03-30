# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("MCMC BN Grzeg Sampling (Fast Tests)")

test_that("2-node Bayesian Network", {
  # set.seed(7101)
  # x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  # x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  # x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  # theData <- data.frame(x1 = x1, x2 = x2,  x3 = x3)
  # 
  # fam <- enumerateBNSpace(3)
  # scores <- logScoreMultDir(fam, theData)
  # 
  # priors <- rep(1/25, 25)
  # scores <- scores - max(scores)
  # expected <- exp(scores)*priors/sum(exp(scores)*priors)
  # 
  # numberOfBurnIn <- 10000
  # numberOfSamples <- 20000
  # 
  # expectedTable <- data.frame(expected = expected * numberOfSamples)
  # row.names(expectedTable) <- lapply(fam, function(network) paste(network, sep = "", collapse = ","))
  # 
  # empty <- list(c(),c(),c())
  # 
  # priorFlat <- function(network) {
  #   1/length(fam)
  # }
  # 
  # sampler <- BNSamplerGrzeg(theData, bn(integer(0), integer(0), integer(0)), priorFlat)
  # samples <- lapply(seq_len(numberOfBurnIn), sampler)
  # samples <- lapply(seq_len(numberOfSamples), sampler)
  # 
  # outTable <- table(factor(unlist(lapply(samples,function(l){
  #   paste(l,sep = "",collapse = ",")}))))
  # 
  # expect_that(as.vector(outTable["integer(0),1,2"]),
  #             is_within(2463, 80))
  # expect_that(as.vector(outTable["2:3,integer(0),integer(0)"]),
  #             is_within(55, 30))
  # expect_that(as.vector(outTable["integer(0),3,1"]),
  #             is_within(879, 35))

})
