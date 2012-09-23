# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick


context("MCMC Gibbs Sampling (Slow Tests)")

# test_that("4-node simulated", {
#   set.seed(5141)
#   x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
#   x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
#   x3 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
#   x4 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
#   theData <- data.frame(x1 = x1, x2 = x2, x3 = x3, x4 = x4)
#   family <- enumerateBNSpace(4)
#   exactScores <- logScoreMultDir(family, theData)
#
#   expectedProbs <- exp(exactScores - logsumexp(exactScores))
#
#   initial <- empty(ncol(theData), "bn")
#   sampler <- BNGibbsSampler(data             = theData,
#                            initial          = initial,
#                            prior            = priorUniform(initial),
#                            maxNumberParents = 3,
#                            moveprobs = c(0, 0, 0, 1))
#
#   samples <- draw(sampler = sampler,
#                  n       = 600,
#                  burnin  = 100,
#                  verbose = T)
#   gibbs <- bnpostmcmc(sampler, samples)
#   actual <- pltabulate(samples)
#
#   exact <- bnpost(family, exactScores, theData)
#   ep(exact)
#
#   expected <- expectedProbs * sum(actual)
#   names(expected) <- as.character(family)
#   o <- match(names(expected), names(actual))
#   actual <- actual[o]
#   actual[is.na(actual)] <- 0
#   d <- cbind(expected = round(expected), actual = actual)
#   d <- transform(d, diff = round(actual - expected))
#
#   chisq.test(actual, p = expectedProbs)
# })
