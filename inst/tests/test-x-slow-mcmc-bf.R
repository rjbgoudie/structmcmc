# Part of the "structural" package, http://github.com/rjbgoudie/structural
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rjbgoudie/structural
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("MCMC BF BN Sampling (Slow Tests)")

test_that("3-node Bayesian Network", {
  set.seed(7101)
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  theData <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  fam <- enumerateBNSpace(3)
  scores <- logScoreMultDir(fam, theData)

  priors <- rep(1/25, 25)
  scores <- scores - max(scores)
  expected <- exp(scores)*priors/sum(exp(scores)*priors)

  numberOfBurnIn <- 10000
  numberOfSamples <- 20000

  expectedTable <- data.frame(expected = expected * numberOfSamples)
  row.names(expectedTable) <- lapply(fam, function(network) paste(network, sep = "", collapse = ","))

  empty <- list(c(),c(),c())

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, bn(integer(0), integer(0), integer(0)), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l)paste(l,sep = "",collapse = ",")))))

  expect_that(as.vector(outTable["2,integer(0),2"]), is_within(2460, 150))
  expect_that(as.vector(outTable["2,integer(0),1"]), is_within(120, 20))
  expect_that(as.vector(outTable["integer(0),1,1"]), is_within(120, 20))
  expect_that(as.vector(outTable["integer(0),c(1,3),1"]), is_within(860, 100))
  expect_that(as.vector(outTable["integer(0),3,1"]), is_within(860, 100))
  expect_that(as.vector(outTable["integer(0),3,integer(0)"]), is_within(1200, 100))
  expect_that(as.vector(outTable["2:3,integer(0),integer(0)"]), is_within(60, 25))
})

test_that("2-node Bayesian Network", {
  set.seed(5432)
  x1 <- as.factor(c(1, 0, 0, 1, 1))
  x2 <- as.factor(c(0, 1, 1, 0, 0))
  theData <- data.frame(x1 = x1, x2 = x2)
  numberOfNodes <- 2

  fam <- enumerateBNSpace(numberOfNodes)
  scores <- logScoreMultDir(fam, theData)

  priors <- rep(1/length(fam), length(fam))
  expected <- exp(scores)*priors/sum(exp(scores)*priors)

  numberOfBurnIn <- 20000
  numberOfSamples <- 20000

  expectedTable <- data.frame(expected = expected * numberOfSamples)
  row.names(expectedTable) <- lapply(fam, function(network) paste(network, sep = "", collapse = ","))

  empty <- list(c(),c(),c())

  priorFlat <- function(network) {
    1/length(fam)
  }

  sampler <- BNSamplerBigFlips(theData, do.call("bn", lapply(seq_len(numberOfNodes),function(i) integer(0))), priorFlat)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l)paste(l,sep = "",collapse = ",")))))

  count_NULL_1 <- as.vector(outTable["integer(0),1"])
  count_2_NULL <- as.vector(outTable["2,integer(0)"])
  count_NULL_NULL <- as.vector(outTable["integer(0),integer(0)"])

  expect_that(count_NULL_1, is_within(9336, 100))
  expect_that(count_2_NULL, is_within(9336, 100))
  expect_that(count_NULL_NULL, is_within(1328, 70))
})

test_that("Nick Podd's example", {

  LDAGUnique <- function(LDAG){
    N <- length(LDAG)
    UDAG <- lapply(1:N,function(x)integer(0))
    for(i in 1:N)
    {
      UDAG[[i]] <- sort(LDAG[[i]])
    }
    UDAG
  }

  TableDAGCheck <- function(LDAG,LTable){
    N <- length(LDAG)
    M <- length(LTable)
    Counter <- 0
    for(i in 1:N)
    {
      Counter <- Counter + (length(LTable[[i]]) == 2^length(LDAG[[i]]))
    }
    Error <- c((Counter == N) , (N==M) , checkAcyclic(LDAG))
    if(!Error[1])stop("Table != DAG")
    if(!Error[2])stop("Table Length != DAG Length")
    if(!Error[3])stop("Not Acyclic")
    !!(prod(Error))
  }

  tmp.orderFunc <- function(tmp.S,tmp.A){
    tmp.len <- length(tmp.S)
    tmp.tabsize <- colSums(tmp.A)
    if (sum(tmp.tabsize)!=0)
    {
      tmp.O <- (1:tmp.len)[tmp.tabsize==0]
      tmp.Order <- (tmp.S)[tmp.tabsize==0]
      tmp.newA <- as.matrix(as.matrix(tmp.A[-tmp.O,])[,-tmp.O])
      tmp.newO <- (tmp.S)[tmp.tabsize!=0]
      tmp.newOrder <- tmp.orderFunc(tmp.newO,tmp.newA)
      tmp.r <- c(tmp.Order,tmp.newOrder)
    }
    else
    {
    tmp.r <- tmp.S
    }
    return(tmp.r)
  }

  BinSampler <- function(LDAG,LTable,N){
    DAGLen <- length(LDAG)
    LDAG <- LDAGUnique(LDAG)
    class(LDAG) <- c("bn", "parental")
    ##Table Check
    if(!TableDAGCheck(LDAG,LTable)) stop("Table/DAG Problem")

    AdjMat <- as.adjacency(LDAG)
    class(LDAG) <- "list"

    States <- 1:length(LDAG)
    Order <- tmp.orderFunc(States,AdjMat)

    SampleTable <- matrix(0,ncol=DAGLen,nrow=N)
    RandomTable <- matrix(runif(N*DAGLen),ncol=DAGLen,nrow=N)

    for( i in (Order))
    {
      LTableLen <- length(LDAG[[i]])
      if(LTableLen==0)
      {
        SampleTable[,i] <- (  RandomTable[,i] < LTable[[i]][1])
      }
      else
      { 

        ts <- SampleTable[,States*AdjMat[,i]]
        tm <- as.matrix(2^(0:(LTableLen-1)))
        multi <- (ts %*% tm) + 1
        rej <- LTable[[i]][multi]
        SampleTable[,i] <- (  RandomTable[,i] < rej)

      }
    }

    return(data.frame(SampleTable))
  }

  DataSize <- 10000
  LDAG <- list(integer(0),integer(0),c(1,2),integer(0))
  LTable <- list(c(0.5), c(0.5), c(0.1,0.5,0.5,0.9), c(0.5)) # Graph parameters
  set.seed(2301)

  ##Generate Data
  theData <- BinSampler(LDAG,LTable,DataSize)
  ## MCMC Initialisation
  Prrior <- function(dag) 1
  MCSize <- 1000
  MSSamples <- 1
  initial <- do.call("bn", lapply(1:length(LDAG),function(i) integer(0)))

  numberOfBurnIn <- 20000
  numberOfSamples <- 20000

  theData <- data.frame(lapply(theData, factor))
  sampler <- sampler <- BNSamplerBigFlips(theData, initial, Prrior)
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  numberOfNodes <- 4
  fam <- enumerateBNSpace(numberOfNodes)
  scores <- logScoreMultDir(fam, theData)

  priors <- 1
  scores <- scores - max(scores)
  expected <- exp(scores)*priors/sum(exp(scores)*priors)

  expectedTable <- data.frame(expected = expected * numberOfSamples)
  row.names(expectedTable) <- lapply(fam, function(network) paste(network, sep = "", collapse = ","))

  outTable <- table(factor(unlist(lapply(samples,function(l)paste(l,sep = "",collapse = ",")))))

  count_NULL_NULL_12_NULL <- as.vector(outTable["integer(0),integer(0),1:2,integer(0)"])
  count_NULL_NULL_12_3 <- as.vector(outTable["integer(0),integer(0),1:2,3"])
  count_NULL_NULL_NULL_NULL <- as.vector(outTable["integer(0),integer(0),integer(0),integer(0)"])
  if (is.na(count_NULL_NULL_NULL_NULL)) count_NULL_NULL_NULL_NULL <- 0

  expect_that(count_NULL_NULL_12_NULL, is_within(17535, 500))
  expect_that(count_NULL_NULL_12_3, is_within(327, 290))
  expect_that(count_NULL_NULL_NULL_NULL, is_within(0, 10))
})
