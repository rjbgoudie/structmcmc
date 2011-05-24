# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("MCMC BN Gibbs Sampling (Fast Tests)")

test_that("Simple test", {
  set.seed(9501)
  dat <- data.frame(x1 = as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0)),
                    x2 = as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0)),
                    x3 = as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0)))

  mcmc <- posterior(data = dat, method = "gibbs", verbose = F,
                    nSamples = 1000, nBurnin = 500)
  exact <- posterior(data = dat, method = "exact", verbose = F)

  epmcmc <- ep(mcmc)
  epexact <- ep(exact)

  expect_that(max(epmcmc - epexact) < 0.05, is_true())
  expect_identical(epmcmc, ep(mcmc, method = "tabulate"))
  expect_identical(epmcmc, ep(mcmc, method = "flatten"))
})

test_that("2-node Bayesian Network", {
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
  numberOfSamples <- 50000

  expectedTable <- data.frame(expected = expected * numberOfSamples)
  row.names(expectedTable) <- lapply(fam, function(network){
    paste(network, sep = "", collapse = ",")
  })

  empty <- list(c(),c(),c())

  priorFlat <- function(network) {
    1/length(fam)
  }

  initial <- bn(integer(0), integer(0), integer(0))
  sampler <- BNGibbsSampler(data             = theData,
                            initial          = initial,
                            prior            = priorFlat,
                            maxNumberParents = 2)

  samples <- draw(sampler,
                  numberOfSamples,
                  burnin = numberOfBurnIn,
                  verbose = F)
  outTable <- table(factor(unlist(lapply(samples,function(l){
    paste(l,sep = "",collapse = ",")}))))

  # expect_that(as.vector(outTable["integer(0),1,2"]),
  #             is_within(2463, 40))
  # expect_that(as.vector(outTable["2:3,integer(0),integer(0)"]),
  #             is_within(55, 30))
  # expect_that(as.vector(outTable["integer(0),3,1"]),
  #             is_within(879, 35))

})

test_that("4-node example", {
  set.seed(7101)
  dat <- data.frame(x1 = sample.int(3, size = 100, replace = T),
                    x2 = sample.int(2, size = 100, replace = T),
                    x3 = sample.int(2, size = 100, replace = T),
                    x4 = sample.int(2, size = 100, replace = T))
  dat <- intAsFDF(dat)
  fam <- enumerateBNSpace(4)
  scores <- logScoreMultDir(fam, dat)
  epost <- bnpost(bnspace = fam, logScore = scores, data = dat)

  numberOfBurnIn <- 10000
  numberOfSamples <- 20000
  initial <- bn(integer(0), integer(0), integer(0), integer(0))
  sampler <- BNGibbsSampler(data             = dat,
                       initial          = initial,
                       prior            = function(x) 1,
                       maxNumberParents = 3)

  samples <- draw(sampler,
                  numberOfSamples,
                  burnin = numberOfBurnIn,
                  verbose = F)
  mpost <- bnpostmcmc(sampler, samples)

  ep(epost)
  ep(mpost)
})

test_that("10-node Bayesian Network", {
  # set.seed(7101)
  # 
  # net <- bn(integer(0), 4L, integer(0), c(22L, 40L), 41L,
  #     22L, integer(0), integer(0), 32L, integer(0), 8L, c(9L,
  #     13L), 1L, c(6L), c(29L, 28L), c(43L, 34L), c(13L,
  #     22L), c(4L, 14L), c(6L, 8L), integer(0), 38L, integer(0),
  #     c(13L, 25L), c(22L, 26L), c(4L, 45L), 14L, c(14L,
  #     31L), integer(0), integer(0), 26L, c(4L, 37L), c(6L, 40L,
  #     45L), integer(0), 45L, c(1L, 37L), c(2L, 7L), c(4L,
  #     40L), 13L, integer(0), integer(0), c(38L, 10L), integer(0), c(1L,
  #     40L), c(1L, 31L), 4L)
  # 
  # root <- as.table(array(c(
  #       0.5, # prob of 1
  #       0.5  # prob of 2
  #     ), 2))
  # 
  # oneParent <- as.table(array(c(
  #             # prob of 1 then 2 given
  #   0.8, 0.2, # p = 1
  #   0.2, 0.8  # p = 2
  # ), c(2, 2)))
  # 
  # twoParents <- as.table(array(c(
  #             # prob of 1 then 2 given
  #   0.8, 0.2, # p1 = 1, p2 = 1
  #   0.2, 0.8, # p1 = 2, p2 = 1
  #   0.2, 0.8, # p1 = 1, p2 = 2
  #   0.2, 0.8  # p1 = 2, p2 = 2
  # ), c(2, 2, 2)))
  # 
  # threeParents <- as.table(array(c(
  #             # prob of 1 then 2 given
  #   0.8, 0.2, # p1 = 1, p2 = 1, p3 = 1
  #   0.2, 0.8, # p1 = 2, p2 = 1, p3 = 1
  #   0.2, 0.8, # p1 = 1, p2 = 2, p3 = 1
  #   0.2, 0.8, # p1 = 2, p2 = 2, p3 = 1
  #   0.2, 0.8, # p1 = 1, p2 = 1, p3 = 2
  #   0.2, 0.8, # p1 = 2, p2 = 1, p3 = 2
  #   0.2, 0.8, # p1 = 1, p2 = 2, p3 = 2
  #   0.2, 0.8  # p1 = 2, p2 = 2, p3 = 2
  # ), c(2, 2, 2, 2)))
  # 
  # cpts <- list(root, oneParent, twoParents, threeParents)
  # 
  # numberOfParents <- sapply(net, length)
  # 
  # cpt <- lapply(seq_along(net), function(i){
  #   cpts[[numberOfParents[i] + 1]]
  # })
  # 
  # theData <- simulate(net, nsim = 300, seed = 1234, ptables = cpt)
  # 
  # priorFlat <- function(network) {
  #   1
  # }
  # numberOfBurnIn <- 100
  # initial <- sampleBN(45, 2)
  # sampler <- BNGibbsSampler(theData,
  #                           initial,
  #                           priorFlat,
  #                           maxNumberParents = 4)
  # samples <- lapply(seq_len(numberOfBurnIn), function(i){
  #   if (i %% 1 == 0){
  #     cat(i, ", ")
  #   }
  #   sampler(i)
  # })
  # 
  # numberOfSamples <- 100
  # 
  # progress <- txtProgressBar(max = numberOfSamples, style = 3)
  # setTxtProgressBar(progress, 0)
  # prog <- 0
  # samples <- lapply(seq_len(numberOfSamples), function(i){
  #   prog <<- prog + 1
  #   setTxtProgressBar(progress, prog)
  #   sampler(i)
  # })
  # outTable <- table(factor(unlist(lapply(samples,function(l){
  #   paste(l,sep = "",collapse = ",")}))))
  #  
  #  
  # expect_that(as.vector(outTable["integer(0),1,2"]),
  #             is_within(2463, 40))
  # expect_that(as.vector(outTable["2:3,integer(0),integer(0)"]),
  #             is_within(55, 30))
  # expect_that(as.vector(outTable["integer(0),3,1"]),
  #             is_within(879, 35))

})
