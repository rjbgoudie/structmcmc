# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

# context("MCMC BN PT Sampling (Fast Tests)")
#
#
#
# test_that("2-node Bayesian Network", {
#   set.seed(7101)
#   x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
#   x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
#   x3 <- as.factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
#   theData <- data.frame(x1 = x1, x2 = x2,  x3 = x3)
#  
#   fam <- enumerateBNSpace(3)
#   scores <- logScoreMultDir(fam, theData)
#  
#   priors <- rep(1/25, 25)
#   scores <- scores - max(scores)
#   expected <- exp(scores)*priors/sum(exp(scores)*priors)
#  
#   numberOfBurnIn <- 1000
#   numberOfSamples <- 50000
#  
#   expectedTable <- data.frame(expected = expected * numberOfSamples)
#   row.names(expectedTable) <- lapply(fam, function(network) paste(network, sep = "", collapse = ","))
#  
#   empty <- list(c(),c(),c())
#
#   priorFlat <- function(network) {
#     1/length(fam)
#   }
#  
#   sampler <- BNSamplerIndep(theData, bn(integer(0), integer(0), integer(0)), priorFlat, maxNumberParents = 2)
#   samples <- lapply(seq_len(numberOfBurnIn), function(i){
#     if (i %% 1000 == 0){
#       cat(i, ", ")
#     }
#     sampler(i)
#   })
#   samples <- lapply(seq_len(numberOfSamples), function(i){
#     if (i %% 1000 == 0){
#       cat(i, ", ")
#     }
#     sampler(i)
#   }) 
#   outTable <- table(factor(unlist(lapply(samples,function(l){
#     paste(l,sep = "",collapse = ",")}))))
#    
#    
#     sampler0 <- BNSampler(theData, empty(3, class = "bn"), priorFlat)
#     samples0 <- draw(sampler0, numberOfSamples)
#
#     mpost0 <- bnpostmcmc(sampler0, samples0)
#     mpostep <- ep(mpost0)
#  
#   expect_that(as.vector(outTable["integer(0),1,2"]),
#               is_within(2463, 40))
#   expect_that(as.vector(outTable["2:3,integer(0),integer(0)"]),
#               is_within(55, 30))
#   expect_that(as.vector(outTable["integer(0),3,1"]),
#               is_within(879, 35))
#
# })
#
# test_that("2-node Bayesian Network", {
#   set.seed(7101)
#   x1 <- as.factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
#   x2 <- as.factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
#   x3 <- as.factor(sample(c(0, 1), size = 100, replace = T))
#   x4 <- as.factor(sample(c(0, 1), size = 100, replace = T))
#   theData <- data.frame(x1 = x1, x2 = x2,  x3 = x3, x4 = x4)
#  
#   fam <- enumerateBNSpace(4)
#   scores <- logScoreMultDir(fam, theData)
#  
#   priors <- rep(1/25, 25)
#   scores <- scores - max(scores)
#   expected <- exp(scores)*priors/sum(exp(scores)*priors)
#  
#   numberOfBurnIn <- 10000
#   numberOfSamples <- 50000
#  
#   expectedTable <- data.frame(expected = expected * numberOfSamples)
#   row.names(expectedTable) <- lapply(fam, function(network) paste(network, sep = "", collapse = ","))
#  
#   empty <- list(c(),c(),c())
#
#   priorFlat <- function(network) {
#     1/length(fam)
#   }
#  
#  
#   sampler <- BNSamplerIndep(theData, bn(integer(0), integer(0), integer(0), integer(0)), priorFlat, maxNumberParents = NULL)
#   samples <- lapply(seq_len(numberOfBurnIn), function(i){
#     if (i %% 1000 == 0){
#       cat(i, ", ")
#     }
#     sampler(i)
#   })
#   samples <- lapply(seq_len(numberOfSamples), function(i){
#     if (i %% 1000 == 0){
#       cat(i, ", ")
#     }
#     sampler(i)
#   }) 
#   outTable <- table(factor(unlist(lapply(samples,function(l){
#     paste(l,sep = "",collapse = ",")}))))
#  
#   expect_that(as.vector(outTable["integer(0),1,2"]),
#               is_within(2463, 40))
#   expect_that(as.vector(outTable["2:3,integer(0),integer(0)"]),
#               is_within(55, 30))
#   expect_that(as.vector(outTable["integer(0),3,1"]),
#               is_within(879, 35))
#
# })
#
#
# test_that("Maybe independent", {
#   # go inside a sampler first
#  
#   m[,] <- 1
#   propL <- vector("list", length = 10000)
#   for (propLi in 1:10000){
#     emptyNetwork <- currentNetwork
#     emptyNetwork[[1]] <- empty(numberOfNodes, class = "bn")
#     emptyNetwork[[2]] <- routes(emptyNetwork[[1]])
#     emptyNetwork[[4]] <- as.adjacency(emptyNetwork[[1]])
#     proposalNetwork <- emptyNetwork
#    
#    
#     nodesConstrainedEmpty <- which(colSums(constraint) == 1 - numberOfNodes)
#     nodesDone <- nodesConstrainedEmpty
#    
#     isCyclic <- F
#    
#     # m is a matrix
#     mt <- m
#     logprob <- 0
#    
#     while (length(nodesDone) != numberOfNodes){
#       # loop over each node
#    
#       ppsu <- which(!is.na(mt), arr.ind = T)
#       mtNotNA <- mt[!is.na(mt)]
#       partition <- logsumexp(1 * mtNotNA)
#       logweights <- 1 * mtNotNA - partition
#       weights <- exp(logweights)
#    
#       w <- sample.int(nrow(ppsu),
#                       size = 1,
#                       prob = weights)
#       node <- ppsu[w, "row"]
#       newparents <- sort.int(possibleParentsList[[node]][[ppsu[w, "col"]]])
#    
#       # update proposalNetwork to be M_{\crosscircle
#    
#       if (sum(proposalNetwork[[2]][node, newparents]) > 0){
#         isCyclic <- T
#         proposalNetwork <- emptyNetwork
#         nodesDone <- nodesConstrainedEmpty
#         logprob <- 0
#         mt <- m
#       }
#       else {
#      #   cat("Plus", logweights[w], "\n")
#         logprob <- logprob + logweights[w]
#         nodesDone <- c(nodesDone, node)
#         proposalNetwork[[1]][[node]] <- newparents
#         proposalNetwork[[4]][newparents, node] <- 1
#         for (i in newparents){
#           proposalNetwork[[2]] <- routesAddEdge(proposalNetwork[[2]], i, node)
#         }
#         mt[node, ] <- rep(NA, maxPossibleParentsList)
#       }
#     }
#     if (propLi %% 100 == 0) cat(propLi, ",")
#     propL[[propLi]] <- proposalNetwork[[1]]
#   }
#  
#   outTable <- table(factor(unlist(lapply(propL,function(l){
#     paste(l,sep = "",collapse = ",")}))))
#  
# })
#
# test_that("Maybe independent", {
#   # nine node example
#  
#   net <- bn(integer(0), integer(0), c(1, 6, 9), c(1, 12), integer(0),
#             c(2, 9), integer(0), c(4, 9, 11), integer(0), integer(0),
#             c(10, 12), 6)
#  
#   root <- as.table(array(c(
#         0.5, # prob of 1
#         0.5  # prob of 2
#       ), 2))
#
#   oneParent <- as.table(array(c(
#               # prob of 1 then 2 given
#     0.8, 0.2, # p = 1
#     0.2, 0.8  # p = 2
#   ), c(2, 2)))
#
#   twoParents <- as.table(array(c(
#               # prob of 1 then 2 given
#     0.8, 0.2, # p1 = 1, p2 = 1
#     0.2, 0.8, # p1 = 2, p2 = 1
#     0.2, 0.8, # p1 = 1, p2 = 2
#     0.2, 0.8  # p1 = 2, p2 = 2
#   ), c(2, 2, 2)))
#
#   threeParents <- as.table(array(c(
#               # prob of 1 then 2 given
#     0.8, 0.2, # p1 = 1, p2 = 1, p3 = 1
#     0.2, 0.8, # p1 = 2, p2 = 1, p3 = 1
#     0.2, 0.8, # p1 = 1, p2 = 2, p3 = 1
#     0.2, 0.8, # p1 = 2, p2 = 2, p3 = 1
#     0.2, 0.8, # p1 = 1, p2 = 1, p3 = 2
#     0.2, 0.8, # p1 = 2, p2 = 1, p3 = 2
#     0.2, 0.8, # p1 = 1, p2 = 2, p3 = 2
#     0.2, 0.8  # p1 = 2, p2 = 2, p3 = 2
#   ), c(2, 2, 2, 2)))
#
#   cpts <- list(root, oneParent, twoParents, threeParents)
#
#   numberOfParents <- sapply(net, length)
#
#   cpt <- lapply(seq_along(net), function(i){
#     cpts[[numberOfParents[i] + 1]]
#   })
#
#   theData <- simulate.bn(net, cpt, N = 1000)
#  
#  
#   priorFlat <- function(network) {
#     1
#   }
#  
#   numberOfBurnIn <- 1000
#   numberOfSamples <- 10000
#  
#   sampler <- BNSamplerIndep(theData, empty(12, class = "bn"), priorFlat, maxNumberParents = 2)
#   samples <- lapply(seq_len(numberOfBurnIn), function(i){
#     if (i %% 5 == 0){
#       cat(i, ", ")
#     }
#     sampler(i)
#   })
#   samples <- lapply(seq_len(numberOfSamples), function(i){
#     if (i %% 1000 == 0){
#       cat(i, ", ")
#     }
#     sampler(i)
#   })
#   outTable <- table(factor(unlist(lapply(samples,function(l){
#     paste(l,sep = "",collapse = ",")}))))
#    
#   ep(z[[1]])
#  
#   sampler0 <- BNSampler(theData, empty(12, class = "bn"), priorFlat)
#   samples0 <- draw(sampler0, numberOfSamples)
#  
#   mpost0 <- bnpostmcmc(sampler0, samples0)
#   mpostep <- ep(mpost0)
#  
# system.time(z <- BNImportanceSampler(theData, priorFlat, maxNumberParents = 3, numberOfSamples = 10000))
#
# outTable <- table(factor(unlist(lapply(samples0,function(l){
#   paste(l,sep = "",collapse = ",")}))))
#
# })
