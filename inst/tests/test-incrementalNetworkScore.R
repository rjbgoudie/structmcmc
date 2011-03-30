# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("incrementalNetworkScore()")

test_that("Cooper & Herskovits 1992", {
  ## Test 1
  # Cooper & Herskovits 1992
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  data <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  numberOfNodes <- 3
  nodesSeq <- seq_len(numberOfNodes)

  network <- bn.list(
    bn(numeric(0), 1, 2),
    bn(numeric(0), 1, numeric(0)),
    bn(numeric(0), c(1, 3), numeric(0)),
    bn(numeric(0), numeric(0), numeric(0)),
    bn(numeric(0), numeric(0), c(1, 2)),
    bn(numeric(0), numeric(0), c(2)),
    bn(numeric(0), 1, c(2, 1))
  )

  numberToTest <- length(network)
  testSeq <- seq_len(numberToTest)
  diffSeq <- seq_len(numberToTest - 1)

  offline <- logScoreMultDir(network, data)
  offlinediff <- sapply(diffSeq, function(i) offline[[i + 1]] - offline[[i]])

  # incremental
  data <- data.frame(lapply(data, as.factor))
  data <- data.matrix(data) - 1
  nl <- apply(data, 2, function(i) length(unique(i)))
  names(nl) <- nodesSeq
  cache <- new.env(hash=TRUE, size = 10000L)

  for (head in nodesSeq){
    localLogScoreMultDir(node = head, parents = network[[1]][[head]], logScoreParameters = list(data = data, nl = nl), cache, checkInput = F)
  }

  diff <- vector("numeric", length = 10)
  diff[1] <- logScoreMultDirIncremental(network[[1]], network[[2]], 3, cache, logScoreParameters = list(data = data, nl = nl), checkInput = F)
  diff[2] <- logScoreMultDirIncremental(network[[2]], network[[3]], 2, cache, logScoreParameters = list(data = data, nl = nl), checkInput = F)
  diff[3] <- logScoreMultDirIncremental(network[[3]], network[[4]], 2, cache, logScoreParameters = list(data = data, nl = nl), checkInput = F)
  diff[4] <- logScoreMultDirIncremental(network[[4]], network[[5]], 3, cache, logScoreParameters = list(data = data, nl = nl), checkInput = F)
  diff[5] <- logScoreMultDirIncremental(network[[5]], network[[6]], 3, cache, logScoreParameters = list(data = data, nl = nl), checkInput = F)
  diff[6] <- logScoreMultDirIncremental(network[[6]], network[[7]], c(2,3), cache, logScoreParameters = list(data = data, nl = nl), checkInput = F)

  sapply(diffSeq, function(i) expect_that(offlinediff[[i]], equals(diff[i])))
})

test_that("Cooper & Herskovits 1992, hyp = .9", {
  ## Test 1
  # Cooper & Herskovits 1992
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  data <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  numberOfNodes <- 3
  nodesSeq <- seq_len(numberOfNodes)

  network <- bn.list(
    bn(numeric(0), 1, 2),
    bn(numeric(0), 1, numeric(0)),
    bn(numeric(0), c(1, 3), numeric(0)),
    bn(numeric(0), numeric(0), numeric(0)),
    bn(numeric(0), numeric(0), c(1, 2)),
    bn(numeric(0), numeric(0), c(2)),
    bn(numeric(0), 1, c(2, 1))
  )

  numberToTest <- length(network)
  testSeq <- seq_len(numberToTest)
  diffSeq <- seq_len(numberToTest - 1)

  offline <- logScoreMultDir(network, data, hyperparameters = "point9")
  offlinediff <- sapply(diffSeq, function(i) offline[[i + 1]] - offline[[i]])

  # incremental
  data <- data.frame(lapply(data, as.factor))
  data <- data.matrix(data) - 1
  nl <- apply(data, 2, function(i) length(unique(i)))
  names(nl) <- nodesSeq
  cache <- new.env(hash=TRUE, size = 10000L)

  for (head in nodesSeq){
    localLogScoreMultDir(node = head, parents = network[[1]][[head]], logScoreParameters = list(data = data, nl = nl, hyperparameters = "point9"), cache, checkInput = F)
  }

  diff <- vector("numeric", length = 10)
  diff[1] <- logScoreMultDirIncremental(network[[1]], network[[2]], 3, cache, logScoreParameters = list(data = data, nl = nl, hyperparameters = "point9"), checkInput = F)
  diff[2] <- logScoreMultDirIncremental(network[[2]], network[[3]], 2, cache, logScoreParameters = list(data = data, nl = nl, hyperparameters = "point9"), checkInput = F)
  diff[3] <- logScoreMultDirIncremental(network[[3]], network[[4]], 2, cache, logScoreParameters = list(data = data, nl = nl, hyperparameters = "point9"), checkInput = F)
  diff[4] <- logScoreMultDirIncremental(network[[4]], network[[5]], 3, cache, logScoreParameters = list(data = data, nl = nl, hyperparameters = "point9"), checkInput = F)
  diff[5] <- logScoreMultDirIncremental(network[[5]], network[[6]], 3, cache, logScoreParameters = list(data = data, nl = nl, hyperparameters = "point9"), checkInput = F)
  diff[6] <- logScoreMultDirIncremental(network[[6]], network[[7]], c(2,3), cache, logScoreParameters = list(data = data, nl = nl, hyperparameters = "point9"), checkInput = F)

  sapply(diffSeq, function(i) expect_that(offlinediff[[i]], equals(diff[i])))
})