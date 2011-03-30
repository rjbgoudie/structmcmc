# Part of the "structural" package, http://github.com/rbtgde/structural
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structural
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("MCMC BF BN Sampling (Fast Tests)")

test_that("2-node Bayesian Network", {
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

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")}))))

  expect_that(as.vector(outTable["integer(0),1"]),
              is_within(326, 35))
  expect_that(as.vector(outTable["2,integer(0)"]),
              is_within(326, 35))
  expect_that(as.vector(outTable["integer(0),integer(0)"]),
              is_within(347, 20))
})

test_that("2-node Bayesian Network, all zeros", {
  x1 <- as.factor(c(0, 0, 0, 0))
  x2 <- as.factor(c(0, 0, 0, 0))
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
    paste(network, sep = "", collapse = ",")
  })

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]),
              is_within(333, 20))
  expect_that(as.vector(outTable["2,integer(0)"]),
              is_within(333, 20))
  expect_that(as.vector(outTable["integer(0),integer(0)"]),
              is_within(333, 20))
})

test_that("2-node Bayesian Network, all 8s", {
  x1 <- as.factor(c(8, 8, 8, 8))
  x2 <- as.factor(c(8, 8, 8, 8))
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
    paste(network, sep = "", collapse = ",")
  })

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l, sep = "",collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]),
              is_within(333, 20))
  expect_that(as.vector(outTable["2,integer(0)"]),
              is_within(333, 20))
  expect_that(as.vector(outTable["integer(0),integer(0)"]),
              is_within(333, 20))
})

test_that("2-node Bayesian Network, non contiguous factor levels", {
  x1 <- as.factor(c(0, 8, 8, 8))
  x2 <- as.factor(c(0, 8, 8, 8))
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
    paste(network, sep = "", collapse = ",")
  })

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l, sep = "", collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]),
              is_within(431, 35))
  expect_that(as.vector(outTable["2,integer(0)"]),
              is_within(431, 35))
  expect_that(as.vector(outTable["integer(0),integer(0)"]),
              is_within(138, 35))
})

test_that("2-node Bayesian Network, non contiguous factor levels", {
  x1 <- as.factor(c(0, 1, 3, 4))
  x2 <- as.factor(c(0, 2, 3, 4))
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
    paste(network, sep = "", collapse = ",")
  })

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]), is_within(433, 35))
  expect_that(as.vector(outTable["2,integer(0)"]), is_within(433, 35))
  expect_that(as.vector(outTable["integer(0),integer(0)"]), is_within(132, 35))
})

test_that("2-node Bayesian Network, non contiguous factor levels", {
  x1 <- as.factor(c(-5:5))
  x2 <- as.factor(c(5:-5))
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
    paste(network, sep = "", collapse = ",")
  })

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]), is_within(494, 35))
  expect_that(as.vector(outTable["2,integer(0)"]), is_within(494, 35))
  expect_that(as.vector(outTable["integer(0),integer(0)"]), is_within(10, 35))
})

test_that("2-node Bayesian Network, non contiguous factor levels", {
  x1 <- as.factor(c(-50:50))
  x2 <- as.factor(c(50:-50))
  theData <- data.frame(x1 = x1, x2 = x2)
  set.seed(1538)
  fam <- enumerateBNSpace(2)
  scores <- logScoreMultDir(fam, theData)

  priors <- 1
  scores <- scores - max(scores)
  expected <- exp(scores)*priors/sum(exp(scores)*priors)

  numberOfBurnIn <- 100
  numberOfSamples <- 1000

  expectedTable <- data.frame(expected = expected * numberOfSamples)
  row.names(expectedTable) <- lapply(fam, function(network){
    paste(network, sep = "", collapse = ",")
  })

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]), is_within(500, 35))
  expect_that(as.vector(outTable["2,integer(0)"]), is_within(500, 35))
})

test_that("2-node Bayesian Network, with random noise", {
  set.seed(1020)
  a <- sample(c(-50:50), size = 50, replace = T)
  b <- sample(c(-50:50), size = 50, replace = T)

  x1 <- as.factor(c(-50:50, a))
  x2 <- as.factor(c(50:-50, b))
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
    paste(network, sep = "", collapse = ",")
  })

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]), is_within(481, 35))
  expect_that(as.vector(outTable["2,integer(0)"]), is_within(481, 35))
  expect_that(as.vector(outTable["integer(0),integer(0)"]), is_within(37, 35))
})

test_that("2-node Bayesian Network, with random noise", {
  set.seed(1020)
  a <- sample(c(-50:50), size = 55, replace = T)
  b <- sample(c(-50:50), size = 55, replace = T)

  x1 <- as.factor(c(-50:50, a))
  x2 <- as.factor(c(50:-50, b))
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
    paste(network, sep = "", collapse = ",")
  })

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]), is_within(469, 35))
  expect_that(as.vector(outTable["2,integer(0)"]), is_within(469, 35))
  expect_that(as.vector(outTable["integer(0),integer(0)"]), is_within(61, 35))
})

test_that("2-node Bayesian Network", {
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
    paste(network, sep = "", collapse = ",")
  })

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]),
              is_within(326, 35))
  expect_that(as.vector(outTable["2,integer(0)"]),
              is_within(326, 35))
  expect_that(as.vector(outTable["integer(0),integer(0)"]),
              is_within(347, 20))

  TrueProbs <- function(theData,type=2){
    N<-length(theData)
    graphs <- enumerateBNSpace(N)
    LenG <- length(graphs)
    z<-numeric(0)
    for(i in 1:LenG)
    {
    if(type==1) z[i] <- LogProbProp(graphs[[i]],theData)
    if(type==2) z[i] <- logScoreMultDir(graphs[[i]], theData)[[1]]
    }
    z <- z -max(z)
    return(list(exp(z)/sum(exp(z)),graphs))
  }

  expect_that(TrueProbs(theData)[[1]], equals(expected))
})

test_that("Non-factor input", {
  x1 <- c(1, 1, 0, 1)
  x2 <- c(0, 1, 0, 1)
  theData <- data.frame(x1 = x1, x2 = x2)

  numberOfBurnIn <- 100
  numberOfSamples <- 1000

  priorFlat <- function(network) {
    1/length(fam)
  }

  expect_that(BNSamplerBigFlips(theData, bn(integer(0), integer(0)), priorFlat),
              throws_error())
})
