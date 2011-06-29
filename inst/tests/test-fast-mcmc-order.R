# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("MCMC BN Order Sampling (Fast Tests)")

test_that("Simple test", {
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
  numberOfSamples <- 50000

  expectedTable <- data.frame(expected = expected * numberOfSamples)
  row.names(expectedTable) <- lapply(fam, function(network){
    paste(network, sep = "", collapse = ",")
  })

  empty <- list(c(),c(),c())

  priorFlat <- function(network) {
    1
  }

  initial <- bn(integer(0), integer(0), integer(0))
  sampler <- BNOrderSampler(data             = theData,
                            initial          = initial,
                            prior            = priorFlat,
                            maxNumberParents = 2)

  order <- 1:3
  numberOfNodes <- get("numberOfNodes", env = environment(sampler))
  nodesSeq <- get("nodesSeq", env = environment(sampler))
  scoresParents <- get("scoresParents", env = environment(sampler))
  parentsTables <- get("parentsTables", env = environment(sampler))
  allRows <- get("allRows", env = environment(sampler))
  rowsThatContain <- get("rowsThatContain", env = environment(sampler))

  actual <- logOrderLikelihood(order           = order,
                               numberOfNodes   = numberOfNodes,
                               nodesSeq        = nodesSeq,
                               scoresParents   = scoresParents,
                               parentsTables   = parentsTables,
                               allRows         = allRows,
                               rowsThatContain = rowsThatContain)

  expected <- log(sum(exp(scores[sapply(fam, isConsistentWithOrder,
                                        order = order)])))

  expect_equal(actual, expected)
})


test_that("Simple test", {
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
  numberOfSamples <- 100000

  expectedTable <- data.frame(expected = expected * numberOfSamples)
  row.names(expectedTable) <- lapply(fam, function(network){
    paste(network, sep = "", collapse = ",")
  })

  empty <- list(c(),c(),c())

  priorFlat <- function(network) {
    1
  }

  initial <- bn(integer(0), integer(0), integer(0))
  sampler <- BNOrderSampler(data             = theData,
                            initial          = 1:3,
                            prior            = priorFlat,
                            maxNumberParents = 2)

  orders <- lapply(seq_len(numberOfBurnIn), sampler)
  orders <- lapply(seq_len(numberOfSamples/50), sampler)

  ellisWong(orders, sampler)

  samples <- unlist(lapply(orders, function(order){
    graphs <- list()

    dagGivenOrder(order,
                    numberOfNodes,
                    nodesSeq,
                    scoresParents,
                    parentsTables,
                    allRows,
                    rowsThatContain)
  }), rec = F)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")}))))

  u <- unique(samples)

  invisible(lapply(u, function(x){

    count <- sum(sapply(samples, function(i){
      identical(i, x)
    }))

    adjusted <- 1/numberOrdersGivenDAG(x) * count
    cat(as.character(x), adjusted, "\n")
  }))
  expect_equal(actual, expected)
})

