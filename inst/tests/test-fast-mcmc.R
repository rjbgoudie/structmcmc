# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("MCMC BN Sampling (Fast Tests)")

test_that("Simple test", {
  set.seed(9501)
  dat <- data.frame(x1 = as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0)),
                    x2 = as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0)),
                    x3 = as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0)))

  mcmc <- posterior(data = dat, method = "mc3", verbose = F,
                    nSamples = 1000, nBurnin = 500)
  exact <- posterior(data = dat, method = "exact", verbose = F)

  epmcmc <- ep(mcmc)
  epexact <- ep(exact)

  expect_that(max(epmcmc - epexact) < 0.05, is_true())
  expect_identical(epmcmc, ep(mcmc, method = "tabulate"))
  expect_identical(epmcmc, ep(mcmc, method = "flatten"))
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
    paste(network, sep = "", collapse = ",")})

  sampler <- BNSampler(theData, bn(integer(0), integer(0)), prior = priorUniform(bn(integer(0), integer(0))))
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

  initial <- empty(2, "bn")

  sampler <- BNSampler(theData, initial, priorUniform(initial))
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

  initial <- empty(2, "bn")

  sampler <- BNSampler(theData, initial, priorUniform(initial))
  samples <- lapply(seq_len(numberOfBurnIn), sampler)
  samples <- lapply(seq_len(numberOfSamples), sampler)

  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l, sep = "",collapse = ",")
  }))))

  expect_that(as.vector(outTable["integer(0),1"]),
              is_within(333, 24))
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

  initial <- empty(2, "bn")

  sampler <- BNSampler(theData, initial, priorUniform(initial))
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

  initial <- empty(2, "bn")

  sampler <- BNSampler(theData, initial, priorUniform(initial))
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

  initial <- empty(2, "bn")

  sampler <- BNSampler(theData, initial, priorUniform(initial))
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

  initial <- empty(2, "bn")

  sampler <- BNSampler(theData, initial, priorUniform(initial))
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

  initial <- empty(2, "bn")

  sampler <- BNSampler(theData, initial, priorUniform(initial))
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

  initial <- empty(2, "bn")

  sampler <- BNSampler(theData, initial, priorUniform(initial))
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

  initial <- empty(2, "bn")

  sampler <- BNSampler(theData, initial, priorUniform(initial))
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

  initial <- empty(2, "bn")

  expect_that(BNSampler(theData, initial, priorUniform(initial)),
              throws_error())
})

test_that("MJ", {
  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  dat <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  bnspace <- enumerateBNSpace(3, allowCyclic = TRUE)
  bnspace <- filterCyclic(bnspace)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  post <- bnpost(bnspace, lsmd, dat)

  pgp <- gp(post)

  nSamples <- 850
  prior <- function(net) 1
  initial <- bn(integer(0), integer(0), integer(0))

  # get "[][1][1]"  "[][1,3][]" "[][][]"
  modes <- bnspace[c(15, 11, 1)]
  class(modes) <- c("bn.list", "parental.list")
  sampler1 <- BNSamplerMJ(dat,
                        initial,
                        prior,
                        modejumping = list(modes = modes,
                                           modeJumpingProbability = 0.25),
                        keepTape = TRUE)

  samples1 <- draw(sampler1, nSamples, verbose = F)

  # plus initial sometimes
  ncgism <- sum(as.vector(sapply(samples1, function(sample){
    sapply(modes, function(mode) identical(sample, mode))}))) + 1

  expectedDiagnostics <- structure(list(
                                 nAccepted        = 530,
                                 nSteps           = 850,
                                 acceptanceRate   = 0.623529411764706,
                                 nMHProposals     = 838,
                                 nMJProposals     = 12,
                                 nMHAccepted      = 527,
                                 nMJAccepted      = 3,
                                 MHAcceptanceRate = 0.628878281622912,
                                 MJAcceptanceRate = 0.25,
                                 nCurrentGraphIsAMode = ncgism),
                          .Names = c("nAccepted",
                                     "nSteps",
                                     "acceptanceRate",
                                     "nMHProposals",
                                     "nMJProposals",
                                     "nMHAccepted",
                                     "nMJAccepted",
                                     "MHAcceptanceRate",
                                     "MJAcceptanceRate",
                                     "nCurrentGraphIsAMode"))

  expect_that(sampler1(returnDiagnostics = T),
              equals(expectedDiagnostics))

  tp <- sampler1(returnTape = T)
  tpdf <- as.data.frame(tp)
  tpdf$gr <- unlist(lapply(samples1, as.character, pretty = T))

  summary(pgp)

  expect_that(
    tpdf[1, "logAccProb"],
    equals(-0.300829647719857)
    # = logScoreMultDir(as.bn("[3][][]", pretty = T), data = dat,
    # hyperparameters = "qi")[[1]] - logScoreMultDir(as.bn("[][][]",
    # pretty = T), data = dat, hyperparameters = "qi")[[1]]
  )

  expect_that(
    tpdf[111, "logAccProb"],
    equals(
      -2.71350405178467
      # = logScoreMultDir(as.bn("[][][]", pretty = T), data = dat,
      # hyperparameters = "qi")[[1]] - logScoreMultDir(as.bn("[][1,3][]",
      #pretty = T), data = dat, hyperparameters = "qi")[[1]] + log(6) - log(6)
    )
  )

  # logScoreMultDir(as.bn("[3][][2]", pretty = T), data = dat,
  # hyperparameters = "qi")[[1]] - logScoreMultDir(as.bn("[3][3][]",
  # pretty = T), data = dat, hyperparameters = "qi")[[1]] + log(6) - log(5)
  expect_that(tpdf[70, "logAccProb"],
              equals(0.182321556793955))

  # logScoreMultDir(as.bn("[2][][]", pretty = T), data = dat,
  # hyperparameters = "qi")[[1]] - logScoreMultDir(as.bn("[2][3][]",
  # pretty = T), data = dat, hyperparameters = "qi")[[1]] + log(5) - log(6)
  expect_that(tpdf[843, "logAccProb"],
              equals(-2.92601434679752))

  # logScoreMultDir(as.bn("[][3][]", pretty = T), data = dat,
  # hyperparameters = "qi")[[1]] - logScoreMultDir(as.bn("[2][3][]",
  # pretty = T), data = dat, hyperparameters = "qi")[[1]] + log(5) - log(6)
  expect_that(tpdf[844, "logAccProb"],
              equals(-0.911111326255255))

  # logScoreMultDir(as.bn("[][1][1]", pretty = T), data = dat,
  # hyperparameters = "qi")[[1]] - logScoreMultDir(as.bn("[][1,3][]",
  # pretty = T), data = dat, hyperparameters = "qi")[[1]] + log(6) - log(6)
  expect_that(tpdf[692, "logAccProb"],
              equals(-2.28554393004323))

})

test_that("Constraints basics", {
  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  dat <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  nSamples <- 850
  # Does not move out of constraint 1
  initial <- bn(integer(0), integer(0), integer(0))
  constraint <- matrix(c( 0, -1, -1,
                         -1,  0, -1,
                          0, -1,  0), nrow = 3, ncol = 3, byrow = T)
  sampler1 <- BNSampler(dat,
                        initial,
                        priorUniform(initial),
                        constraint = constraint)

  samples1 <- draw(sampler1, nSamples, verbose = F)

  allInAllowedSet <- all(sapply(samples1, function(x){
    isOption1 <- identical(x, bn(3L, integer(0), integer(0)))
    isOption2 <- identical(x, bn(integer(0), integer(0), integer(0)))
    isOption1 | isOption2
  }))
  expect_that(allInAllowedSet, is_identical_to(TRUE))

  # Does not move out of constraint 2
  initial <- bn(3L, integer(0), integer(0))
  constraint <- matrix(c( 0, -1, -1,
                         -1,  0, -1,
                          1,  0,  0), nrow = 3, ncol = 3, byrow = T)
  sampler1 <- BNSampler(dat,
                        initial,
                        priorUniform(initial),
                        constraint = constraint)

  samples1 <- draw(sampler1, nSamples, verbose = F)

  allInAllowedSet <- all(sapply(samples1, function(x){
    isOption1 <- identical(x, bn(3L, integer(0), integer(0)))
    isOption2 <- identical(x, bn(3L, 3L, integer(0)))
    isOption1 | isOption2
  }))
  expect_that(allInAllowedSet, is_identical_to(TRUE))

  # Does not move out of constraint 3
  initial <- bn(integer(0), c(1L), c(1L, 2L))
  constraint <- matrix(c( 0,  1,  0,
                         -1,  0,  1,
                         -1, -1,  0), nrow = 3, ncol = 3, byrow = T)
  sampler1 <- BNSampler(dat,
                        initial,
                        priorUniform(initial),
                        constraint = constraint)

  samples1 <- draw(sampler1, nSamples, verbose = F)

  allInAllowedSet <- all(sapply(samples1, function(x){
    isOption1 <- identical(x, bn(integer(0), c(1L), c(1L, 2L)))
    isOption2 <- identical(x, bn(integer(0), c(1L), c(2L)))
    isOption1 | isOption2
  }))
  expect_that(allInAllowedSet, is_identical_to(TRUE))

  # Initial does not meet constraint 1
  initial <- bn(integer(0), c(1L), integer(0))
  constraint <- matrix(c( 0, -1, -1,
                         -1,  0, -1,
                         -1, -1,  0), nrow = 3, ncol = 3, byrow = T)

  expect_that(BNSampler(dat,
                        initial,
                        priorUniform(initial),
                        constraint = constraint),
              throws_error("Initial network does not satisfy constraint"))

  # Initial does not meet constraint 2
  initial <- bn(integer(0), 1L, c(2L))
  constraint <- matrix(c( 0,  1,  1,
                          0,  0,  1,
                          0,  0,  0), nrow = 3, ncol = 3, byrow = T)

  expect_that(BNSampler(dat,
                        initial,
                        priorUniform(initial),
                        constraint = constraint),
              throws_error("Initial network does not satisfy constraint"))

  # Initial does not meet constraint 3
  initial <- bn(integer(0), 1L, c(2L))
  constraint <- matrix(c( 0, -1,  1,
                          0,  0,  1,
                          0,  0,  0), nrow = 3, ncol = 3, byrow = T)

  expect_that(BNSampler(dat,
                        initial,
                        priorUniform(initial),
                        constraint = constraint),
              throws_error("Initial network does not satisfy constraint"))

  # Gets correct probabilities 1
  set.seed(7201)
  constraint <- matrix(c( 0, -1, -1,
                         -1,  0, -1,
                          0, -1,  0), nrow = 3, ncol = 3, byrow = T)

  bnspace <- enumerateBNSpace(3)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  data.frame(lsmd, as.character(bnspace))

  # allowed graphs are
  # 23 == [3][][]
  # 25 == [][][]
  data.frame(lsmd, as.character(bnspace))[c(23, 25), ]

  # so probs should be, given flat prior
  normalisingConstant <- logsumexp(lsmd[c(23, 25)])
  expectedProbs <- exp(lsmd[c(23, 25)] - normalisingConstant)
  expectedNums <- expectedProbs * nSamples

  initial <- bn(integer(0), integer(0), integer(0))
  sampler1 <- BNSampler(dat,
                        initial,
                        priorUniform(initial),
                        constraint = constraint)

  samples1 <- draw(sampler1, nSamples, verbose = F)

  outTable <- table(factor(unlist(lapply(samples1, function(l){
    paste(l, sep = "", collapse = ",")}))))

  expect_that(as.vector(outTable["3,integer(0),integer(0)"]),
              is_within(361, 10))

  expect_that(as.vector(outTable["integer(0),integer(0),integer(0)"]),
              is_within(488, 10))

  # Gets correct probabilities 2
  constraint <- matrix(c( 0, -1, -1,
                         -1,  0, -1,
                          1,  0,  0), nrow = 3, ncol = 3, byrow = T)

  set.seed(4301)
  bnspace <- enumerateBNSpace(3)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  data.frame(lsmd, as.character(bnspace))

  # allowed graphs are
  # 17 == [3][3][]
  # 23 == [3][][]
  data.frame(lsmd, as.character(bnspace))[c(17, 23), ]

  # so probs should be, given flat prior
  normalisingConstant <- logsumexp(lsmd[c(17, 23)])
  expectedProbs <- exp(lsmd[c(17, 23)] - normalisingConstant)
  expectedNums <- expectedProbs * nSamples

  initial <- bn(3L, integer(0), integer(0))
  sampler1 <- BNSampler(dat,
                        initial,
                        priorUniform(initial),
                        constraint = constraint)

  samples1 <- draw(sampler1, nSamples, verbose = F)

  outTable <- table(factor(unlist(lapply(samples1, function(l){
    paste(l, sep = "", collapse = ",")}))))

  expect_that(as.vector(outTable["3,3,integer(0)"]),
              is_within(799, 10))

  expect_that(as.vector(outTable["3,integer(0),integer(0)"]),
              is_within(51, 10))
})

test_that("Using Prior for indegree constraint", {
  set.seed(9501)
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  dat <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  nSamples <- 850
  
  localPriors <- lapply(seq_len(ncol(dat)), function(i){
    function(x){
      if (length(x) > 1){
        0
      } else {
        1
      }
    }
  })

  bnspace <- enumerateBNSpace(3)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  data.frame(lsmd, as.character(bnspace))

  whichTwoParents <- c(3, 9, 11, 12, 13, 18, 20, 21, 24)
  # allowed graphs are
  data.frame(lsmd, as.character(bnspace))[-whichTwoParents, ]

  # so probs should be, given flat prior
  normalisingConstant <- logsumexp(lsmd[-whichTwoParents])
  expectedProbs <- exp(lsmd[-whichTwoParents] - normalisingConstant)
  expectedNums <- expectedProbs * nSamples
  data.frame(gr = as.character(bnspace)[-whichTwoParents], expectedNums)

  initial <- bn(integer(0), integer(0), integer(0))
  sampler1 <- BNSampler(dat,
                        initial,
                        prior = localPriors)
  samples1 <- draw(sampler1, nSamples, verbose = F)

  outTable <- table(factor(unlist(lapply(samples1, function(l){
    paste(l, sep = "", collapse = ",")}))))

  expect_that(as.vector(outTable["2,3,integer(0)"]),
              is_within(157, 15))

  expect_that(as.vector(outTable["integer(0),1,2"]),
              is_within(157, 15))

 expect_that(as.vector(outTable["integer(0),1,1"]),
             is_within(8, 10))
})

test_that("Proper indegree constraint", {
  set.seed(9501)
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  dat <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  nSamples <- 850

  bnspace <- enumerateBNSpace(3)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  data.frame(lsmd, as.character(bnspace))

  whichTwoParents <- c(3, 9, 11, 12, 13, 18, 20, 21, 24)
  # allowed graphs are
  data.frame(lsmd, as.character(bnspace))[-whichTwoParents, ]

  # so probs should be, given flat prior
  normalisingConstant <- logsumexp(lsmd[-whichTwoParents])
  expectedProbs <- exp(lsmd[-whichTwoParents] - normalisingConstant)
  expectedNums <- expectedProbs * nSamples
  data.frame(gr = as.character(bnspace)[-whichTwoParents], expectedNums)

  initial <- bn(integer(0), integer(0), integer(0))
  sampler1 <- BNSampler(dat,
                        initial,
                        prior = priorUniform(initial),
                        maxNumberParents = 1)
  samples1 <- draw(sampler1, nSamples, verbose = F)
  samples1 <- draw(sampler1, nSamples, verbose = F)

  outTable <- table(factor(unlist(lapply(samples1, function(l){
    paste(l, sep = "", collapse = ",")}))))

  expect_that(as.vector(outTable["2,3,integer(0)"]),
              is_within(157, 20))

  expect_that(as.vector(outTable["integer(0),1,2"]),
              is_within(157, 20))

 expect_that(as.vector(outTable["integer(0),1,1"]),
             is_within(8, 10))
})

test_that("Mode-jumping with constraint", {
  set.seed(310)
  x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
  x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
  x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
  dat <- data.frame(x1 = x1, x2 = x2,  x3 = x3)

  nSamples <- 850
  initial <- bn(3L, integer(0), integer(0))
  constraint <- matrix(c( 0,  0, -1,
                         -1,  0, -1,
                          1,  0,  0), nrow = 3, ncol = 3, byrow = T)

  bnspace <- enumerateBNSpace(3)
  lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
  data.frame(lsmd, as.character(bnspace))

  # allowed graphs are
  # 6 == [3][1][]
  # 9 == [3][3][]
  # 12 == [3][1,3][]
  # 3 == [3][][]
  data.frame(lsmd, as.character(bnspace))[c(6, 9, 12, 3), ]

  # so probs should be, given flat prior
  normalisingConstant <- logsumexp(lsmd[c(6, 9, 12, 3)])
  expectedProbs <- exp(lsmd[c(6, 9, 12, 3)] - normalisingConstant)
  expectedNums <- expectedProbs * nSamples
  data.frame(gr = as.character(bnspace)[c(6, 9, 12, 3)], expectedNums)

  modes <- bnspace[c(12, 3)]
  class(modes) <- c("bn.list", "parental.list")
  sampler1 <- BNSamplerMJ(dat,
                        initial,
                        prior = priorUniform(initial),
                        modejumping = list(modes = modes,
                                           modeJumpingProbability = 0.25),
                        constraint = constraint)

  samples1 <- draw(sampler1, nSamples, verbose = F)

  outTable <- table(factor(unlist(lapply(samples1, function(l){
    paste(l, sep = "", collapse = ",")}))))

  expect_that(as.vector(outTable["3,1,integer(0)"]),
              is_within(52, 16))
  expect_that(as.vector(outTable["3,c(1,3),integer(0)"]),
              is_within(380, 25))
  expect_that(as.vector(outTable["3,integer(0),integer(0)"]),
              is_within(25, 16))
})
