# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

# context("SHD Plot")
#
# test_that("Very weird input", {
#  
#   set.seed(5391)
#  
#   computeScores <- function(N){
#     net <- bn(integer(0), c(1,3), 4, integer(0))
#
#     root <- as.table(array(c(
#     0.5, # prob of 1
#     0.5  # prob of 2
#     ), 2))
#
#     oneParent <- as.table(array(c(
#     # prob of 1 then 2 given
#     0.8, 0.2, # p = 1
#     0.2, 0.8  # p = 2
#     ), c(2, 2)))
#
#     twoParents <- as.table(array(c(
#                 # prob of 1 then 2 given
#       0.8, 0.2, # p1 = 1, p2 = 1
#       0.2, 0.8, # p1 = 2, p2 = 1
#       0.2, 0.8, # p1 = 1, p2 = 2
#       0.2, 0.8  # p1 = 2, p2 = 2
#     ), c(2, 2, 2)))
#
#     cpts <- list(root, oneParent, twoParents)
#
#     numberOfParents <- sapply(net, length)
#
#     cpt <- lapply(seq_along(net), function(i){
#       cpts[[numberOfParents[i] + 1]]
#     })
#
#     dat <- simulate.bn(net, cpt, N)
#    
#     fam <- enumerateBNSpace(length(net))
#     scores <- logScoreMultDir(fam, dat)
#     epost <- bnpost(fam, scores, dat)
#     tp <- top(epost)
#  
#     tpScores <- logScoreMultDir(tp, dat)
#  
#     tpSeq <- seq_along(tp)
#     tpSize <- length(tp)
#     m <- matrix(NA, nrow = tpSize, ncol = tpSize)
#  
#     for (i in tpSeq){
#       for (j in tpSeq){
#         m[i, j] <- numberOfMovesBetweenIgnoringCycles(x = tp[[i]],
#                                                       y = tp[[j]],
#                                                       allowFlips = T)
#       }
#     }
#     colnames(m) <- as.character(tp, pretty = T)
#     rownames(m) <- as.character(tp, pretty = T)
#  
#     c(sum(m), sum(abs(1 - outer(tpScores, 1/tpScores))))
#   }
#  
#  
#   Ns <- c(300, 500, 1000, 2000)
#   sapply(Ns, computeScores)
#  
#   levelplot(m,
#             scales = list(x = list(rot = 90)))
# })
