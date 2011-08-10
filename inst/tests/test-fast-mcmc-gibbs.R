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
                    nSamples = 50000, nBurnin = 10000)
  mcmcmh <- posterior(data = dat, method = "mc3", verbose = F,
                    nSamples = 50000, nBurnin = 10000)
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
  theData <- data.frame(x1 = x1, x2 = x2)
  family <- enumerateBNSpace(2)
  exactScores <- logScoreMultDir(family, theData)

  expectedProbs <- exp(scores - logsumexp(exactScores))

  sampler <- BNGibbsSampler(data             = theData,
                            initial          = empty(ncol(theData), "bn"),
                            prior            = priorUniform(),
                            maxNumberParents = 2,
                            moveprobs = c(0, 1, 0))

  samples <- draw(sampler = sampler,
                  n       = 50000,
                  burnin  = 10000,
                  verbose = T)
  gibbs <- bnpostmcmc(sampler, samples)
  actual <- pltabulate(samples)

  expected <- expectedProbs * sum(actual)
  names(expected) <- as.character(family)
  o <- match(names(expected), names(actual))
  actual <- actual[o]
  actual[is.na(actual)] <- 0
  d <- cbind(expected = round(expected), actual = actual)
  d <- transform(d, diff = round(actual - expected))
  
  chisq.test(actual, p = expectedProbs)



  # expect_that(as.vector(outTable["integer(0),1,2"]),
  #             is_within(2463, 40))
  # expect_that(as.vector(outTable["2:3,integer(0),integer(0)"]),
  #             is_within(55, 30))
  # expect_that(as.vector(outTable["integer(0),3,1"]),
  #             is_within(879, 35))

})



test_that("5-node Bayesian Network", {
  #set.seed(7101)
  set.seed(5141)
  dat <- data.frame(x1 = sample.int(3, size = 100, replace = T),
                    x2 = sample.int(2, size = 100, replace = T),
                    x3 = sample.int(2, size = 100, replace = T),
                    x4 = sample.int(2, size = 100, replace = T),
                    x5 = sample.int(2, size = 100, replace = T))
  theData <- intAsFDF(dat)
  family <- enumerateBNSpace(ncol(theData))
  exactScores <- logScoreMultDir(family, theData)

  expectedProbs <- exp(exactScores - logsumexp(exactScores))


  expected <- expectedProbs * sum(actual)
  names(expected) <- as.character(family)
  o <- match(names(expected), names(actual))
  actual <- actual[o]
  actual[is.na(actual)] <- 0
  d <- cbind(expected = round(expected), actual = actual)
  d <- transform(d, diff = round(actual - expected))
  
  chisq.test(actual, p = expectedProbs)

  sampler2 <- BNSampler(data             = theData,
                            initial          = empty(ncol(theData), "bn"),
                            prior            = priorUniform(),
                            maxNumberParents = 5)

  samples2 <- draw(sampler = sampler2,
                  n       = 50000,
                  burnin  = 10000,
                  verbose = T)
  mh <- bnpostmcmc(sampler2, samples2)
  
  
  
  sampler <- BNGibbsSampler(data             = theData,
                            initial          = empty(ncol(theData), "bn"),
                            prior            = priorUniform(),
                            maxNumberParents = 5,
                            moveprobs = c(0, 1, 0))

  samples <- draw(sampler = sampler,
                  n       = 1000,
                  burnin  = 10,
                  verbose = T)
  actual <- pltabulate(samples)
  gibbs <- bnpostmcmc(sampler, samples)
  exact <- bnpost(family, exactScores, theData)


  samples <- draw(sampler = gsampler,
                  n       = 1000,
                  burnin  = 10,
                  verbose = T)

  sampler3 <- BNGibbsSampler(data             = theData,
                             initial          = empty(ncol(theData), "bn"),
                             prior            = priorUniform(),
                             maxNumberParents = 5,
                             moveprobs = c(0, 0.9, 0.1))

   samples3 <- draw(sampler = sampler3,
                   n       = 50000,
                   burnin  = 10000,
                   verbose = T)
   gibbs3 <- bnpostmcmc(sampler3, samples3)
   
   sampler4 <- BNGibbsSampler(data             = theData,
                              initial          = empty(ncol(theData), "bn"),
                              prior            = priorUniform(),
                              maxNumberParents = 5,
                              moveprobs = c(0.75, 0.2, 0.05))

    samples4 <- draw(sampler = sampler4,
                    n       = 50000,
                    burnin  = 10000,
                    verbose = T)
    gibbs4 <- bnpostmcmc(sampler4, samples4)

    sampler5 <- BNGibbsSampler(data             = theData,
                               initial          = empty(ncol(theData), "bn"),
                               prior            = priorUniform(),
                               maxNumberParents = 5,
                               moveprobs = c(0, 0, 1))

     samples5 <- draw(sampler = sampler5,
                     n       = 50000,
                     burnin  = 10000,
                     verbose = T)
     gibbs5 <- bnpostmcmc(sampler5, samples5)

     xyplot(cumtvd(gp(exact), bnpostmcmc.list(gibbs3 = gibbs3, gibbs5 = gibbs5)))

  xyplot(cumtvd(gp(exact), bnpostmcmc.list(gibbs = gibbs, mh = mh, gibbs3=gibbs3, gibbs4=gibbs4)))


  # expect_that(as.vector(outTable["integer(0),1,2"]),
  #             is_within(2463, 40))
  # expect_that(as.vector(outTable["2:3,integer(0),integer(0)"]),
  #             is_within(55, 30))
  # expect_that(as.vector(outTable["integer(0),3,1"]),
  #             is_within(879, 35))

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
  expectedProbs <- exp(scores)*priors/sum(exp(scores)*priors)

  sampler <- BNGibbsSampler(data             = theData,
                            initial          = empty(ncol(dat), "bn"),
                            prior            = priorUniform(),
                            maxNumberParents = 2,
                            moveprobs = c(0, 0, 1))

  samples <- draw(sampler = sampler,
                  n       = 5000,
                  burnin  = 1000,
                  verbose = T)
  actual <- pltabulate(samples)

  expected <- expectedProbs * sum(actual)
  names(expected) <- as.character(fam)
  o <- match(names(expected), names(actual))
  actual <- actual[o]
  actual[is.na(actual)] <- 0
  d <- cbind(expected = round(expected), actual = actual)
  d <- transform(d, diff = round(actual - expected))

  chisq.test(actual, p = expectedProbs)


  # expect_that(as.vector(outTable["integer(0),1,2"]),
  #             is_within(2463, 40))
  # expect_that(as.vector(outTable["2:3,integer(0),integer(0)"]),
  #             is_within(55, 30))
  # expect_that(as.vector(outTable["integer(0),3,1"]),
  #             is_within(879, 35))

})
# 
# test_that("4-node example", {
#   set.seed(7101)
#   dat <- data.frame(x1 = sample.int(3, size = 100, replace = T),
#                     x2 = sample.int(2, size = 100, replace = T),
#                     x3 = sample.int(2, size = 100, replace = T),
#                     x4 = sample.int(2, size = 100, replace = T))
#   dat <- intAsFDF(dat)
#   fam <- enumerateBNSpace(4)
#   scores <- logScoreMultDir(fam, dat)
#   epost <- bnpost(bnspace = fam, logScore = scores, data = dat)
# 
#   numberOfBurnIn <- 10
#   numberOfSamples <- 75
#   initial <- bn(integer(0), integer(0), integer(0), integer(0))
#   sampler <- BNGibbsSampler(data             = dat,
#                        initial          = initial,
#                        prior            = function(x) 1,
#                        maxNumberParents = 3,
#                        moveprobs = c(0, 1, 0),
#                        fam =fam,
#                        exactScore = scores)
# 
#   samples <- draw(sampler,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = T)
#   mpost <- bnpostmcmc(sampler, samples)
#   ep(epost)
#   ep(mpost)
#   
#   
#   numberOfBurnIn <- 100
#   numberOfSamples <- 75000
#   initial <- bn(integer(0), integer(0), integer(0), integer(0))
#   sampler2 <- BNGibbsSampler(data             = dat,
#                        initial          = initial,
#                        prior            = function(x) 1,
#                        maxNumberParents = 3, moveprobs = c(1, 0, 0))
# 
#   samples2 <- draw(sampler2,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = T)
#   mpost2 <- bnpostmcmc(sampler2, samples2)
#   
#   
#   
#   numberOfBurnIn <- 100
#   numberOfSamples <- 1000
#   initial <- bn(integer(0), integer(0), integer(0), integer(0))
#   sampler3 <- BNGibbsSampler(data             = dat,
#                        initial          = initial,
#                        prior            = function(x) 1,
#                        maxNumberParents = 3, moveprobs = c(0, 0, 1))
# 
#   samples3 <- draw(sampler3,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = T)
#   mpost3 <- bnpostmcmc(sampler3, samples3)
#   
#   xyplot(cumtvd(gp(epost), bnpostmcmc.list(mpost, mpost2, mpost3)))
#   
# })
# 
# 
# test_that("5 node 1", {
#   set.seed(711)
#   dat <- data.frame(x1 = sample.int(3, size = 100, replace = T),
#                     x2 = sample.int(2, size = 100, replace = T),
#                     x3 = sample.int(2, size = 100, replace = T),
#                     x4 = sample.int(2, size = 100, replace = T),
#                     x5 = sample.int(4, size = 100, replace = T))
#   dat <- intAsFDF(dat)
#   fam <- enumerateBNSpace(5)
#   scores <- logScoreMultDir(fam, dat)
#   epost <- bnpost(bnspace = fam, logScore = scores, data = dat)
# 
#   numberOfBurnIn <- 100
#   numberOfSamples <- 100000
#   initial <- bn(integer(0), integer(0), integer(0), integer(0), integer(0))
#   sampler <- BNGibbsSampler(data             = dat,
#                        initial          = initial,
#                        prior            = function(x) 1,
#                        maxNumberParents = 4, moveprobs = c(0, 1, 0),
#                        fam = fam,
#                        exactScore = scores)
# 
#   samples <- draw(sampler,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = T)
#   mpost <- bnpostmcmc(sampler, samples)
#   ep(epost)
#   ep(mpost)
#   
#   
#   numberOfBurnIn <- 100
#   numberOfSamples <- 100000
#   initial <- bn(integer(0), integer(0), integer(0), integer(0), integer(0))
#   sampler2 <- BNGibbsSampler(data             = dat,
#                        initial          = initial,
#                        prior            = function(x) 1,
#                        maxNumberParents = 3, moveprobs = c(1, 0, 0))
# 
#   samples2 <- draw(sampler2,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = T)
#   mpost2 <- bnpostmcmc(sampler2, samples2)
#   
#   
#   
#   numberOfBurnIn <- 100
#   numberOfSamples <- 2000
#   initial <- bn(integer(0), integer(0), integer(0), integer(0), integer(0))
#   sampler3 <- BNGibbsSampler(data             = dat,
#                        initial          = initial,
#                        prior            = function(x) 1,
#                        maxNumberParents = 3, moveprobs = c(0, 0, 1))
# 
#   samples3 <- draw(sampler3,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = T)
#   mpost3 <- bnpostmcmc(sampler3, samples3)
#   
#   xyplot(cumtvd(gp(epost), bnpostmcmc.list(mpost, mpost2, mpost3)))
#   
# })
# 
# test_that("3 node 1", {
#   set.seed(711)
#   dat <- data.frame(x1 = sample.int(3, size = 100, replace = T),
#                     x2 = sample.int(2, size = 100, replace = T),
#                     x3 = sample.int(2, size = 100, replace = T))
#   dat <- intAsFDF(dat)
#   fam <- enumerateBNSpace(3)
#   scores <- logScoreMultDir(fam, dat)
#   epost <- bnpost(bnspace = fam, logScore = scores, data = dat)
# 
#   numberOfBurnIn <- 100
#   numberOfSamples <- 100000
#   initial <- bn(integer(0), integer(0), integer(0))
#   sampler <- BNGibbsSampler(data             = dat,
#                        initial          = initial,
#                        prior            = function(x) 1,
#                        maxNumberParents = 3, moveprobs = c(0, 1, 0))
# 
#   samples <- draw(sampler,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = T)
#   mpost <- bnpostmcmc(sampler, samples)
#   ep(epost)
#   ep(mpost)
#   
#   
#   numberOfBurnIn <- 100
#   numberOfSamples <- 100000
#   initial <- bn(integer(0), integer(0), integer(0))
#   sampler2 <- BNGibbsSampler(data             = dat,
#                        initial          = initial,
#                        prior            = function(x) 1,
#                        maxNumberParents = 3, moveprobs = c(1, 0, 0))
# 
#   samples2 <- draw(sampler2,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = T)
#   mpost2 <- bnpostmcmc(sampler2, samples2)
#   
#   
#   
#   numberOfBurnIn <- 100
#   numberOfSamples <- 2000
#   initial <- bn(integer(0), integer(0), integer(0))
#   sampler3 <- BNGibbsSampler(data             = dat,
#                        initial          = initial,
#                        prior            = function(x) 1,
#                        maxNumberParents = 3, moveprobs = c(0, 0, 1))
# 
#   samples3 <- draw(sampler3,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = T)
#   mpost3 <- bnpostmcmc(sampler3, samples3)
#   
#   xyplot(cumtvd(gp(epost), bnpostmcmc.list(mpost, mpost2, mpost3)))
#   
# })
# 
# test_that("5-node example", {
#   source("~/Desktop/fam5.R")
#   library(structmcmc)
#   set.seed(7101)
#   dat <- data.frame(x1 = sample.int(3, size = 100, replace = T),
#                     x2 = sample.int(2, size = 100, replace = T),
#                     x3 = sample.int(2, size = 100, replace = T),
#                     x4 = sample.int(2, size = 100, replace = T),
#                     x5 = sample.int(2, size = 100, replace = T))
#   dat <- intAsFDF(dat)
# 
#   numberOfBurnIn <- 100
#   numberOfSamples <- 2000
#   initial <- empty(5, "bn")
#   sampler <- BNGibbsSampler(data             = dat,
#                             initial          = initial,
#                             prior            = function(x) 1,
#                             maxNumberParents = 4,
#                             moveprobs = c(0, 1, 0))
#   sapply(1:10, sampler)
# 
#   samples <- draw(sampler,
#                   numberOfSamples,
#                   burnin = numberOfBurnIn,
#                   verbose = F)
#   mpost <- bnpostmcmc(sampler, samples)
#   
#   fam <- enumerateBNSpace(5)
#   scores <- logScoreMultDir(fam, dat)
#   epost <- bnpost(bnspace = fam, logScore = scores, data = dat)
#   
#   rows2graph <- function(rows, change, currentNetwork, parentsTables){
#     net <- currentNetwork[[1]]
#     numberOfNodes <- nNodes(net)
#     base <- empty(numberOfNodes, "bn")
#     nodesSeq <- seq_len(numberOfNodes)
#     notchange <- setdiff(nodesSeq, change)
#     base[notchange] <- net[notchange]
#     
#     out <- lapply(rows, function(row){
#       ops <- expand.grid(row)
#       out <- list()
#       for (i in seq_len(nrow(ops))){
#         gr <- base
#         this <- as.matrix(ops[i, ])
#         for (j in seq_along(change)){
#           node <- change[j]
#           new <- parentsTables[[node]][this[j], ]
#           new <- new[!is.na(new)]
#           if (length(new) > 0){
#             gr[[node]] <- new
#           }
#         }
#         out <- c(out, list(gr))
#       }
#       out
#     })
#     out
#   }
#   out <- rows2graph(rows, nodes, currentNetwork, parentsTables)
#   
#   out <- unlist(out, rec = F)
#   
#   net <- currentNetwork[[1]]
#   numberOfNodes <- nNodes(net)
#   nodesSeq <- seq_len(numberOfNodes)
#   notchange <- setdiff(nodesSeq, nodes)
# 
#   fam <- enumerateBNSpace(5)
#   
#   filter <- function(x){
#     if (identical(x[[2]], integer(0)) && identical(x[[4]], integer(0))){
#       T
#     } else {
#       F
#     }
#   }
#   
#   
#   actual <- fam[sapply(fam, filter)]
#   expected <- unlist(out, rec = F)
#   
#   find <- function(x, y){
#     any(sapply(y, function(z){
#       identical(x, z)
#     }))
#   }
#   
#   stopifnot(all(sapply(actual, function(x) find(x, expected))),
#             all(sapply(expected, function(x) find(x, actual))))
# 
#   filter2 <- function(x, net, notchange){
#     is.ok <- T
#     for (i in notchange){
#       if (!identical(net[[i]], x[[i]])){
#         is.ok <- F
#       }
#     }
#     is.ok
#   }
#   
#   sapply(fam, filter2, currentNetwork[[1]], notchange)
# 
#   ep(epost)
#   ep(mpost)
# })



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



  # nodes <- c(node1, node2)
  # rows <- list(list(rows1[[1]], rows2[[1]]),
  #              list(rows1[[2]], rows2[[2]]),
  #              list(rows1[[3]], rows2[[3]]))
  # rows2graph <- function(rows, change, currentNetwork, parentsTables){
  #   net <- currentNetwork[[1]]
  #   numberOfNodes <- nNodes(net)
  #   base <- empty(numberOfNodes, "bn")
  #   nodesSeq <- seq_len(numberOfNodes)
  #   notchange <- setdiff(nodesSeq, change)
  #    base[notchange] <- net[notchange]
  # 
  #   out <- lapply(rows, function(row){
  #     ops <- expand.grid(row)
  #     out <- list()
  #     for (i in seq_len(nrow(ops))){
  #       gr <- base
  #       this <- as.matrix(ops[i, ])
  #       for (j in seq_along(change)){
  #         node <- change[j]
  #         new <- parentsTables[[node]][this[j], ]
  #         new <- new[!is.na(new)]
  #         if (length(new) > 0){
  #           gr[[node]] <- new
  #         }
  #       }
  #       out <- c(out, list(gr))
  #     }
  #     out
  #   })
  #   out
  # }
  # out <- rows2graph(rows, nodes, currentNetwork, parentsTables)
  # 
  # filter2 <- function(x, net, notchange){
  #   is.ok <- T
  #   for (i in notchange){
  #     if (!identical(net[[i]], x[[i]])){
  #       is.ok <- F
  #     }
  #   }
  #   is.ok
  # }
  # 
  # net <- currentNetwork[[1]]
  # numberOfNodes <- nNodes(net)
  # nodesSeq <- seq_len(numberOfNodes)
  # notchange <- setdiff(nodesSeq, nodes)
  # 
  # whFam <- sapply(fam, filter2, currentNetwork[[1]], notchange)
  # expected <- fam[whFam]
  # expectedScores <- exactScore[whFam]
  # expectedScoresN <- exp(expectedScores - logsumexp(expectedScores))
  # actual <- unlist(out, rec = F)
  # 
  # find <- function(x, y){
  #   any(sapply(y, function(z){
  #     identical(x, z)
  #   }))
  # }
  # 
  # if (!isTRUE(all(sapply(expected, function(x) find(x, actual))))){
  #   browser()
  # }
  # 
  # if (!isTRUE(all(sapply(actual, function(x) find(x, expected))))){
  #   browser()
  # }
  # 
  # if (!isTRUE(length(expected) == length(actual))){
  #   browser()
  # }
  # 
  # 
  # sampleGraph <- function(groupWeights, rows1, rows2){
  #   n2SampGroup <- sample.int(3, size = 1, prob = groupWeights)
  # 
  #   # sample 'node1' parents
  #   n1scoresGroup <- scoresParents[[node1]][rows1[[n2SampGroup]]]
  #   n1probs <- exp(n1scoresGroup - logsumexp(n1scoresGroup))
  #   n1samp <- sample.int(length(n1scoresGroup), size = 1, prob = n1probs)
  # 
  #   # sample 'node2' parents
  #   n2scoresGroup <- scoresParents[[node2]][rows2[[n2SampGroup]]]
  #   n2probs <- exp(n2scoresGroup - logsumexp(n2scoresGroup))
  #   n2samp <- sample.int(length(n2scoresGroup), size = 1, prob = n2probs)
  # 
  #   # generate the new graph
  #   parents1 <- rows1[[n2SampGroup]]
  #   new <- parentsTables[[node1]][parents1[n1samp], ]
  #   currentNetwork[[1]][[node1]] <- new[!is.na(new)]
  # 
  #   parents2 <- rows2[[n2SampGroup]]
  #   new <- parentsTables[[node2]][parents2[n2samp], ]
  #   currentNetwork[[1]][[node2]] <- new[!is.na(new)]
  #   currentNetwork[[1]]
  # }
  # 
  # actual <- lapply(1:1000, function(i){
  #   sampleGraph(groupWeights, rows1, rows2)
  # })
  # class(actual) <- "parental.list"
  # expectedNames <- as.character(expected)
  # actualTable <- pltabulate(actual)
  # o <- match(expectedNames, names(actualTable))
  # p <- function(a, p) choose(1000, a) * p ^ a * (1 - p)^(1000 - a)
  # a <- actualTable[o]
  # d <- data.frame(e = expectedScoresN * 1000, a = a)
  # d <- transform(d, diff = round(e - a))
  # a[is.na(a)] <- 0
  # p <- chisq.test(a, p = expectedScoresN)$p.value
  # if (p < 0.05){
  #   browser()
  # }
  #
  
  
  
  # rows2graph <- function(rows, change, currentNetwork, parentsTables){
  #   net <- currentNetwork[[1]]
  #   numberOfNodes <- nNodes(net)
  #   base <- empty(numberOfNodes, "bn")
  #   nodesSeq <- seq_len(numberOfNodes)
  #   notchange <- setdiff(nodesSeq, change)
  #   base[notchange] <- net[notchange]
  #   
  #   out <- lapply(rows, function(row){
  #     ops <- expand.grid(row)
  #     out <- list()
  #     for (i in seq_len(nrow(ops))){
  #       gr <- base
  #       this <- as.matrix(ops[i, ])
  #       for (j in seq_along(change)){
  #         node <- change[j]
  #         new <- parentsTables[[node]][this[j], ]
  #         new <- new[!is.na(new)]
  #         if (length(new) > 0){
  #           gr[[node]] <- new
  #         }
  #       }
  #       out <- c(out, list(gr))
  #     }
  #     out
  #   })
  #   out
  # }
  # out <- rows2graph(rows, nodes, currentNetwork, parentsTables)
  # 
  # filter2 <- function(x, net, notchange){
  #   is.ok <- T
  #   for (i in notchange){
  #     if (!identical(net[[i]], x[[i]])){
  #       is.ok <- F
  #     }
  #   }
  #   is.ok
  # }
  # 
  # net <- currentNetwork[[1]]
  # numberOfNodes <- nNodes(net)
  # nodesSeq <- seq_len(numberOfNodes)
  # notchange <- setdiff(nodesSeq, nodes)
  # 
  # expected <- fam[sapply(fam, filter2, currentNetwork[[1]], notchange)]
  # actual <- unlist(out, rec = F)
  # 
  # find <- function(x, y){
  #   any(sapply(y, function(z){
  #     identical(x, z)
  #   }))
  # }
  # 
  # if (!isTRUE(all(sapply(expected, function(x) find(x, actual))))){
  #   browser()
  # }
  # 
  # if (!isTRUE(all(sapply(actual, function(x) find(x, expected))))){
  #   browser()
  # }
  # 
  # if (!isTRUE(length(expected) == length(actual))){
  #   browser()
  # }
  # 
  # 
  # sampleGraph <- function(groupWeights, rows1, rows2){
  #   sampGroup <- sample.int(25, size = 1, prob = groupWeights)
  #   
  #   # sample 'node1' parents
  #   n1scoresGroup <- scoresParents[[node1]][rows[[sampGroup]][[1]]]
  #   n1probs <- exp(n1scoresGroup - logsumexp(n1scoresGroup))
  #   n1samp <- sample.int(length(n1scoresGroup), size = 1, prob = n1probs)
  # 
  #   # sample 'node2' parents
  #   n2scoresGroup <- scoresParents[[node2]][rows[[sampGroup]][[2]]]
  #   n2probs <- exp(n2scoresGroup - logsumexp(n2scoresGroup))
  #   n2samp <- sample.int(length(n2scoresGroup), size = 1, prob = n2probs)
  # 
  #   # sample 'node3' parents
  #   n3scoresGroup <- scoresParents[[node3]][rows[[sampGroup]][[3]]]
  #   n3probs <- exp(n3scoresGroup - logsumexp(n3scoresGroup))
  #   n3samp <- sample.int(length(n3scoresGroup), size = 1, prob = n3probs)
  # 
  #   # generate the new graph
  #   parents1 <- rows[[sampGroup]][[1]]
  #   new <- parentsTables[[node1]][parents1[n1samp], ]
  #   currentNetwork[[1]][[node1]] <- new[!is.na(new)]
  # 
  #   parents2 <- rows[[sampGroup]][[2]]
  #   new <- parentsTables[[node2]][parents2[n2samp], ]
  #   currentNetwork[[1]][[node2]] <- new[!is.na(new)]
  # 
  #   parents3 <- rows[[sampGroup]][[3]]
  #   new <- parentsTables[[node3]][parents3[n3samp], ]
  #   currentNetwork[[1]][[node3]] <- new[!is.na(new)]
  #   currentNetwork[[1]]
  # }
  # 
  # actual <- lapply(1:1000, function(i){
  #   sampleGraph(groupWeights, rows1, rows2)
  # })
  # class(actual) <- "parental.list"
  # expectedNames <- as.character(expected)
  # actualTable <- pltabulate(actual)
  # o <- match(expectedNames, names(actualTable))
  # p <- function(a, p) choose(1000, a) * p ^ a * (1 - p)^(1000 - a)
  # a <- actualTable[o]
  # d <- data.frame(e = expectedScoresN * 1000, a = a)
  # d <- transform(d, diff = round(e - a))
  # a[is.na(a)] <- 0
  # p <- chisq.test(a, p = expectedScoresN)$p.value
  # if (p < 0.05){
  #   browser()
  # }
  # 
  


})




