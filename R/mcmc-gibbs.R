# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Sample the parents of a single node (Gibbs sampler).
#'
#' Sample from posterior distribution on graph, conditional on
#' all the edges, except those that go into node \code{node}.
#'
#' @param currentNetwork A \code{currentNetwork} object
#' @param numberOfNodes The number of nodes in the network. A numeric vector 
#'   of length 1.
#' @param nodesSeq The vector 1:nNodes(currentNetwork). (Supplied as an 
#'   argument for possible speed gain)
#' @param scoresParents A list of the form returned by 
#'   \code{scoreParentsTable()}
#' @param parentsTables A list of tables of the form returned by 
#'   \code{enumerateParentsTable()}
#' @param allRows The vector 1:nrow(parentsTables). (Supplied as an 
#'   argument for possible speed gain)
#' @param rowsThatContain A list of the form created by 
#'   \code{getRowsThatContain()}
#' @return Returns the sampled network. A \code{currentNetwork} object.
#' @export
#' @seealso \code{\link{BNGibbsSampler}}, \code{\link{samplePair}}
sampleNode <- function(currentNetwork,
                       numberOfNodes,
                       nodesSeq,
                       scoresParents,
                       parentsTables,
                       allRows,
                       rowsThatContain){
  # choose a node to sample the parents of
  node <- sample.int(numberOfNodes, size = 1)

  # remove the old parents of node 'node'
  currentNetwork[[2]] <- routesRemoveEdges(currentNetwork[[2]],
                                           currentNetwork[[1]][[node]],
                                           node)

  # get the conditional probability for the parents of
  # node 'node', given the rest of the graph
  nonDescendants <- nonDescendants(currentNetwork[[2]], node)
  rows <- whichParentSetRows(node            = node,
                             nonDescendants  = nonDescendants,
                             numberOfNodes   = numberOfNodes,
                             allRows         = allRows,
                             rowsThatContain = rowsThatContain)
  scores <- scoresParents[[node]][rows]
  scoresNormalised <- exp(scores - logsumexp(scores))

  # sample a new parent set, according to the condtional probability
  samp <- sample.int(length(scores),
                     size = 1,
                     prob = scoresNormalised)

  # set the new graph to the sampled graph
  new <- parentsTables[[node]][rows[samp], ]
  currentNetwork[[1]][[node]] <- new[!is.na(new)]
  currentNetwork[[2]] <- routesAddEdges(currentNetwork[[2]],
                                        currentNetwork[[1]][[node]],
                                        node)
  currentNetwork[[4]][currentNetwork[[1]][[node]], node] <- 1
  if (length(currentNetwork[[1]][[node]]) == 0){
    currentNetwork[[4]][, node] <- 0
  } else {
    currentNetwork[[4]][-currentNetwork[[1]][[node]], node] <- 0
  }

  currentNetwork[[5]][node] <- scores[samp]

  currentNetwork
}

#' Sample the parents of a pair of nodes (Gibbs sampler).
#'
#' Sample from posterior distribution on graph, conditional on
#' all the edges, except for those corresponding to the parents sets of 
#' two nodes.
#'
#' @param currentNetwork A \code{currentNetwork} object
#' @param numberOfNodes The number of nodes in the network. A numeric vector 
#'   of length 1.
#' @param nodesSeq The vector 1:nNodes(currentNetwork). (Supplied as an 
#'   argument for possible speed gain)
#' @param scoresParents A list of the form returned by 
#'   \code{scoreParentsTable()}
#' @param parentsTables A list of tables of the form returned by 
#'   \code{enumerateParentsTable()}
#' @param allRows The vector 1:nrow(parentsTables). (Supplied as an 
#'   argument for possible speed gain)
#' @param rowsThatContain A list of the form created by 
#'   \code{getRowsThatContain()}
#' @return Returns the sampled network. A \code{currentNetwork} object.
#' @export
#' @seealso \code{\link{BNGibbsSampler}}, \code{\link{sampleNode}}
samplePair <- function(currentNetwork,
                       numberOfNodes,
                       nodesSeq,
                       scoresParents,
                       parentsTables,
                       allRows,
                       rowsThatContain){
  node1 <- sample.int(numberOfNodes, size = 1)
  choices <- setdiff3(nodesSeq, node1)
  node2 <- choices[sample.int(length(choices), size = 1)]

  # remove the old parents of node 'node1' and 'node2'

  currentNetwork[[2]] <- routesRemoveEdges(currentNetwork[[2]],
                                           currentNetwork[[1]][[node1]],
                                           node1)
  currentNetwork[[2]] <- routesRemoveEdges(currentNetwork[[2]],
                                           currentNetwork[[1]][[node2]],
                                           node2)

  nonDescendants1 <- nonDescendants(currentNetwork[[2]], node1)
  descendants1 <- setdiff3(nodesSeq, nonDescendants1)
  nonDescendants2 <- nonDescendants(currentNetwork[[2]], node2)
  descendants2 <- setdiff3(nodesSeq, nonDescendants2)

  rows1 <- vector("list", 3)
  rows2 <- vector("list", 3)

  # group1. nonDescendants1 = nonDescendants1 union nonDescendants2
  newNonDescendants1 <- intersect2(nonDescendants1, nonDescendants2)

  # get the conditional probability for the parents of
  # node 'node1', given the rest of the graph
  rows1[[1]] <- whichParentSetRows(node            = node1,
                                   nonDescendants  = newNonDescendants1,
                                   numberOfNodes   = numberOfNodes,
                                   allRows         = allRows,
                                   rowsThatContain = rowsThatContain)

  # haveNewDescendants == F
  # no new nonDescendants2
  newNonDescendants2 <- intersect2(nonDescendants1, nonDescendants2)

  rows2[[1]] <- whichParentSetRows(node            = node2,
                                   nonDescendants  = newNonDescendants2,
                                   numberOfNodes   = numberOfNodes,
                                   allRows         = allRows,
                                   rowsThatContain = rowsThatContain)
  if (length(rows1[[1]]) > 0 && length(rows2[[1]]) > 0){
    group1Score <- logsumexp(scoresParents[[node1]][rows1[[1]]]) +
                   logsumexp(scoresParents[[node2]][rows2[[1]]])
    group1Score <- sum(group1Score)
  } else {
    group1Score <- -Inf
  }

  # group2
  newNonDescendants1 <- nonDescendants1
  # get the conditional probability for the parents of
  # node 'node1', given the rest of the graph
  rows1[[2]] <- whichParentSetRows(node            = node1,
                                   nonDescendants  = newNonDescendants1,
                                   needOneOf       = descendants2,
                                   numberOfNodes   = numberOfNodes,
                                   allRows         = allRows,
                                   rowsThatContain = rowsThatContain)

  # haveNewDescendants == T
  # nonDescendants2 = nonDescendants2 + descendants1
  # == intersect2(nonDescendants1, nonDescendants2)
  newNonDescendants2 <- intersect2(nonDescendants2, nonDescendants1)
  rows2[[2]] <- whichParentSetRows(node            = node2,
                                   nonDescendants  = newNonDescendants2,
                                   numberOfNodes   = numberOfNodes,
                                   allRows         = allRows,
                                   rowsThatContain = rowsThatContain)
  if (length(rows1[[2]]) > 0 && length(rows2[[2]]) > 0){
    group2Score <- logsumexp(scoresParents[[node1]][rows1[[2]]]) +
                   logsumexp(scoresParents[[node2]][rows2[[2]]])
    group2Score <- sum(group2Score)
  } else {
    group2Score <- -Inf
  }

  # group3
  newNonDescendants1 <- intersect2(nonDescendants1, nonDescendants2)

  # get the conditional probability for the parents of
  # node 'node1', given the rest of the graph
  rows1[[3]] <- whichParentSetRows(node            = node1,
                                   nonDescendants  = newNonDescendants1,
                                   numberOfNodes   = numberOfNodes,
                                   allRows         = allRows,
                                   rowsThatContain = rowsThatContain)

  # haveNewDescendants == T
  # nonDescendants2 = nonDescendants2 + descendants1
  # == intersect2(nonDescendants1, nonDescendants2)
  newNonDescendants2 <- nonDescendants2
  rows2[[3]] <- whichParentSetRows(node            = node2,
                                   nonDescendants  = newNonDescendants2,
                                   needOneOf       = descendants1,
                                   numberOfNodes   = numberOfNodes,
                                   allRows         = allRows,
                                   rowsThatContain = rowsThatContain)
  if (length(rows1[[3]]) > 0 && length(rows2[[3]]) > 0){
    group3Score <- logsumexp(scoresParents[[node1]][rows1[[3]]])
                   logsumexp(scoresParents[[node2]][rows2[[3]]])
    group3Score <- sum(group3Score)
  } else {
    group3Score <- -Inf
  }


  # sample group
  groupScoresOld <- c(group1Score, group2Score, group3Score)
  groupWeights <- exp(groupScoresOld - logsumexp(groupScoresOld))
  n2SampGroup <- sample.int(3, size = 1, prob = groupWeights)

  # sample 'node1' parents
  n1scoresGroup <- scoresParents[[node1]][rows1[[n2SampGroup]]]
  n1probs <- exp(n1scoresGroup - logsumexp(n1scoresGroup))
  n1samp <- sample.int(length(n1scoresGroup), size = 1, prob = n1probs)

  # sample 'node2' parents
  n2scoresGroup <- scoresParents[[node2]][rows2[[n2SampGroup]]]
  n2probs <- exp(n2scoresGroup - logsumexp(n2scoresGroup))
  n2samp <- sample.int(length(n2scoresGroup), size = 1, prob = n2probs)

  # generate the new graph
  parents1 <- rows1[[n2SampGroup]]
  new <- parentsTables[[node1]][parents1[n1samp], ]
  currentNetwork[[1]][[node1]] <- new[!is.na(new)]

  parents2 <- rows2[[n2SampGroup]]
  new <- parentsTables[[node2]][parents2[n2samp], ]
  currentNetwork[[1]][[node2]] <- new[!is.na(new)]

  currentNetwork[[2]] <- routesAddEdges(currentNetwork[[2]],
                                        currentNetwork[[1]][[node1]],
                                        node1)
  currentNetwork[[2]] <- routesAddEdges(currentNetwork[[2]],
                                        currentNetwork[[1]][[node2]],
                                        node2)
  currentNetwork[[4]][currentNetwork[[1]][[node1]], node1] <- 1
  if (length(currentNetwork[[1]][[node1]]) == 0){
    currentNetwork[[4]][, node1] <- 0
  } else {
    currentNetwork[[4]][-currentNetwork[[1]][[node1]], node1] <- 0
  }
  currentNetwork[[4]][currentNetwork[[1]][[node2]], node2] <- 1
  if (length(currentNetwork[[1]][[node2]]) == 0){
    currentNetwork[[4]][, node2] <- 0
  } else {
    currentNetwork[[4]][-currentNetwork[[1]][[node2]], node2] <- 0
  }

  currentNetwork[[5]][node1] <- n1scoresGroup[n1samp]
  currentNetwork[[5]][node2] <- n2scoresGroup[n2samp]

  currentNetwork
}

#' Sample the parents of a pair of nodes (Gibbs sampler) v2.
#'
#' Sample from posterior distribution on graph, conditional on
#' all the edges, except for those corresponding to the parents sets of
#' two nodes.
#'
#' @param currentNetwork A \code{currentNetwork} object
#' @param numberOfNodes The number of nodes in the network. A numeric vector
#'   of length 1.
#' @param nodesSeq The vector 1:nNodes(currentNetwork). (Supplied as an
#'   argument for possible speed gain)
#' @param scoresParents A list of the form returned by
#'   \code{scoreParentsTable()}
#' @param parentsTables A list of tables of the form returned by
#'   \code{enumerateParentsTable()}
#' @param allRows The vector 1:nrow(parentsTables). (Supplied as an
#'   argument for possible speed gain)
#' @param rowsThatContain A list of the form created by
#'   \code{getRowsThatContain()}
#' @return Returns the sampled network. A \code{currentNetwork} object.
#' @export
#' @seealso \code{\link{BNGibbsSampler}}, \code{\link{sampleNode}}
samplePair2 <- function(currentNetwork,
                         numberOfNodes,
                         nodesSeq,
                         scoresParents,
                         parentsTables,
                         allRows,
                         rowsThatContain,
                         logScoreFUN,
                         logScoreParameters){
  node1 <- sample.int(numberOfNodes, size = 1)
  choices <- setdiff3(nodesSeq, node1)
  node2 <- choices[sample.int(length(choices), size = 1)]

  currentNetwork[[2]] <- routesRemoveEdges(currentNetwork[[2]],
                                           currentNetwork[[1]][[node1]],
                                           node1)
  currentNetwork[[2]] <- routesRemoveEdges(currentNetwork[[2]],
                                           currentNetwork[[1]][[node2]],
                                           node2)

  nonDescendants1 <- nonDescendants(currentNetwork[[2]], node1)
  descendants1 <- setdiff3(nodesSeq, nonDescendants1)
  nonDescendants2 <- nonDescendants(currentNetwork[[2]], node2)
  descendants2 <- setdiff3(nodesSeq, nonDescendants2)

  optionsGivenGraph <- function(net){
    if (identical(net[[1]], 2L)){
      newNonDescendants1 <- nonDescendants1
      needOneOf <- descendants2
    } else {
      newNonDescendants1 <- intersect2(nonDescendants1, nonDescendants2)
      needOneOf <- NULL
    }
    rows1 <- whichParentSetRows(node            = node1,
                                nonDescendants  = newNonDescendants1,
                                needOneOf       = needOneOf,
                                numberOfNodes   = numberOfNodes,
                                allRows         = allRows,
                                rowsThatContain = rowsThatContain)

    if (identical(net[[2]], 1L)){
      newNonDescendants2 <- nonDescendants2
      needOneOf <- descendants1
    } else {
      newNonDescendants2 <- intersect2(nonDescendants2, nonDescendants1)
      needOneOf <- NULL
    }
    rows2 <- whichParentSetRows(node            = node2,
                                nonDescendants  = newNonDescendants2,
                                needOneOf       = needOneOf,
                                numberOfNodes   = numberOfNodes,
                                allRows         = allRows,
                                rowsThatContain = rowsThatContain)
    list(rows1, rows2)
  }

  getScoreFromRows <- function(rows){
    if (length(rows[[1]]) > 0 && length(rows[[2]]) > 0){
      groupScore <- logsumexp(scoresParents[[node1]][rows[[1]]])
                    logsumexp(scoresParents[[node2]][rows[[2]]])
      sum(groupScore)
    } else {
      -Inf
    }
  }

  nets <- bn.list(bn(integer(0), integer(0)),
                  bn(2L, integer(0)),
                  bn(integer(0), 1L))
  rows <- lapply(nets, optionsGivenGraph)
  groupScoresOld <- sapply(rows, getScoreFromRows)

  # sample group
  groupWeights <- exp(groupScoresOld - logsumexp(groupScoresOld))
  n2SampGroup <- sample.int(3, size = 1, prob = groupWeights)

  # sample 'node1' parents
  n1scoresGroup <- scoresParents[[node1]][rows[[n2SampGroup]][[1]]]
  n1probs <- exp(n1scoresGroup - logsumexp(n1scoresGroup))
  n1samp <- sample.int(length(n1scoresGroup), size = 1, prob = n1probs)

  # sample 'node2' parents
  n2scoresGroup <- scoresParents[[node2]][rows[[n2SampGroup]][[2]]]
  n2probs <- exp(n2scoresGroup - logsumexp(n2scoresGroup))
  n2samp <- sample.int(length(n2scoresGroup), size = 1, prob = n2probs)

  # generate the new graph
  parents1 <- rows[[n2SampGroup]][[1]]
  new <- parentsTables[[node1]][parents1[n1samp], ]
  currentNetwork[[1]][[node1]] <- new[!is.na(new)]

  parents2 <- rows[[n2SampGroup]][[2]]
  new <- parentsTables[[node2]][parents2[n2samp], ]
  currentNetwork[[1]][[node2]] <- new[!is.na(new)]

  currentNetwork[[2]] <- routesAddEdges(currentNetwork[[2]],
                                        currentNetwork[[1]][[node1]],
                                        node1)
  currentNetwork[[2]] <- routesAddEdges(currentNetwork[[2]],
                                        currentNetwork[[1]][[node2]],
                                        node2)
  currentNetwork[[4]][currentNetwork[[1]][[node1]], node1] <- 1
  if (length(currentNetwork[[1]][[node1]]) == 0){
    currentNetwork[[4]][, node1] <- 0
  } else {
    currentNetwork[[4]][-currentNetwork[[1]][[node1]], node1] <- 0
  }
  currentNetwork[[4]][currentNetwork[[1]][[node2]], node2] <- 1
  if (length(currentNetwork[[1]][[node2]]) == 0){
    currentNetwork[[4]][, node2] <- 0
  } else {
    currentNetwork[[4]][-currentNetwork[[1]][[node2]], node2] <- 0
  }

  currentNetwork[[5]][node1] <- n1scoresGroup[n1samp]
  currentNetwork[[5]][node2] <- n2scoresGroup[n2samp]

  currentNetwork
}


#' Sample the parents of a triple of nodes (Gibbs sampler).
#'
#' Sample from posterior distribution on graph, conditional on
#' all the edges, except for those corresponding to the parents sets of 
#' two nodes.
#'
#' @param currentNetwork A \code{currentNetwork} object
#' @param numberOfNodes The number of nodes in the network. A numeric vector 
#'   of length 1.
#' @param nodesSeq The vector 1:nNodes(currentNetwork). (Supplied as an 
#'   argument for possible speed gain)
#' @param scoresParents A list of the form returned by 
#'   \code{scoreParentsTable()}
#' @param parentsTables A list of tables of the form returned by 
#'   \code{enumerateParentsTable()}
#' @param allRows The vector 1:nrow(parentsTables). (Supplied as an 
#'   argument for possible speed gain)
#' @param rowsThatContain A list of the form created by 
#'   \code{getRowsThatContain()}
#' @param logScoreFUN A list of four elements:
#'   \describe{
#'     \item{offline}{A function that computes the logScore of a Bayesian 
#'                    Network}
#'     \item{online}{A function that incrementally computes the logScore of a 
#'                   Bayesian Network}
#'     \item{local}{A function that computes the local logScore of a 
#'                  Bayesian Network}
#'     \item{prepare}{A function that prepares the data, and any further 
#'                    pre-computation required by the logScore functions.}
#'   }
#'   For Multinomial-Dirichlet models, \code{\link{logScoreMultDirFUN}} 
#'   returns the appropriate list; for Normal models with Zellner g-priors,
#'   \code{\link{logScoreNormalFUN}} returns the appropriate list.
#' @param logScoreParameters A list of parameters that are passed to
#'   logScoreFUN.
#' @return Returns the sampled network. A \code{currentNetwork} object.
#' @export
#' @seealso \code{\link{BNGibbsSampler}}, \code{\link{sampleNode}}
sampleTriple <- function(currentNetwork,
                         numberOfNodes,
                         nodesSeq,
                         scoresParents,
                         parentsTables,
                         allRows,
                         rowsThatContain,
                         logScoreFUN,
                         logScoreParameters){
  node1 <- sample.int(numberOfNodes, size = 1)
  choices <- setdiff3(nodesSeq, node1)
  node2 <- choices[sample.int(length(choices), size = 1)]
  choices <- setdiff3(choices, node2)
  node3 <- choices[sample.int(length(choices), size = 1)]

  nodes <- c(node1, node2, node3)
  currentNetwork[[2]] <- routesRemoveEdges(currentNetwork[[2]],
                                           currentNetwork[[1]][[node1]],
                                           node1)
  currentNetwork[[2]] <- routesRemoveEdges(currentNetwork[[2]],
                                           currentNetwork[[1]][[node2]],
                                           node2)
  currentNetwork[[2]] <- routesRemoveEdges(currentNetwork[[2]],
                                           currentNetwork[[1]][[node3]],
                                           node3)

  # indexed by node1, node2, node3
  # output is for real node numbers
  nonDescendants <- vector("list", 3)
  descendants <- vector("list", 3)
  nonDescendants[[1]] <- nonDescendants(currentNetwork[[2]], node1)
  descendants[[1]] <- setdiff3(nodesSeq, nonDescendants[[1]])
  nonDescendants[[2]] <- nonDescendants(currentNetwork[[2]], node2)
  descendants[[2]] <- setdiff3(nodesSeq, nonDescendants[[2]])
  nonDescendants[[3]] <- nonDescendants(currentNetwork[[2]], node3)
  descendants[[3]] <- setdiff3(nodesSeq, nonDescendants[[3]])

  optionsGivenGraph <- function(net){
    intersectAll <- intersection(nonDescendants[[1]],
                                 nonDescendants[[2]],
                                 nonDescendants[[3]])

    newNonDescendants1 <- intersectAll
    newNonDescendants2 <- intersectAll
    newNonDescendants3 <- intersectAll

    needOneOf1 <- NULL
    needOneOf2 <- NULL
    needOneOf3 <- NULL

    three <- 1:3

    if (length(net[[1]]) > 0){
      descendantOfParent <- descendants[net[[1]]]
      descendantOfNonParent <- unlist(descendants[setdiff3(three, net[[1]])],
                                      use.names = F)
      needOneOf1 <- lapply(descendantOfParent, function(x){
        setdiff3(x, descendantOfNonParent)
      })
      newNonDescendants1 <- c(intersectAll, unlist(needOneOf1, use.names = F))
    }
    if (length(net[[2]]) > 0){
      descendantOfParent <- descendants[net[[2]]]
      descendantOfNonParent <- unlist(descendants[setdiff3(three, net[[2]])],
                                      use.names = F)
      needOneOf2 <- lapply(descendantOfParent, function(x){
        setdiff3(x, descendantOfNonParent)
      })
      newNonDescendants2 <- c(intersectAll, unlist(needOneOf2, use.names = F))
    }
    if (length(net[[3]]) > 0){
      descendantOfParent <- descendants[net[[3]]]
      descendantOfNonParent <- unlist(descendants[setdiff3(three, net[[3]])],
                                      use.names = F)
      needOneOf3 <- lapply(descendantOfParent, function(x){
        setdiff3(x, descendantOfNonParent)
      })
      newNonDescendants3 <- c(intersectAll, unlist(needOneOf3, use.names = F))
    }
    rows1 <- whichParentSetRows(node            = node1,
                                nonDescendants  = newNonDescendants1,
                                needOneOf       = needOneOf1,
                                numberOfNodes   = numberOfNodes,
                                allRows         = allRows,
                                rowsThatContain = rowsThatContain)

    rows2 <- whichParentSetRows(node            = node2,
                                nonDescendants  = newNonDescendants2,
                                needOneOf       = needOneOf2,
                                numberOfNodes   = numberOfNodes,
                                allRows         = allRows,
                                rowsThatContain = rowsThatContain)

    rows3 <- whichParentSetRows(node            = node3,
                                nonDescendants  = newNonDescendants3,
                                needOneOf       = needOneOf3,
                                numberOfNodes   = numberOfNodes,
                                allRows         = allRows,
                                rowsThatContain = rowsThatContain)
    list(rows1, rows2, rows3)
  }

  getScoreFromRows <- function(rows){
    if (length(rows[[1]]) > 0 &&
        length(rows[[2]]) > 0 &&
        length(rows[[3]]) > 0){
      groupScore <- logsumexp(scoresParents[[node1]][rows[[1]]]) +
                    logsumexp(scoresParents[[node2]][rows[[2]]]) +
                    logsumexp(scoresParents[[node3]][rows[[3]]])
      sum(groupScore)
    } else {
      -Inf
    }
  }

  nets <- structure(list(
    structure(list(integer(0), integer(0), integer(0)),
              class = c("bn", "parental")),
    structure(list(2L, integer(0), integer(0)),
              class = c("bn", "parental")),
    structure(list(3L, integer(0), integer(0)),
              class = c("bn", "parental")),
    structure(list(2:3, integer(0), integer(0)),
              class = c("bn","parental")),
    structure(list(integer(0), 1L, integer(0)),
              class = c("bn", "parental")),
    structure(list(3L, 1L, integer(0)),
              class = c("bn", "parental")),
    structure(list(integer(0), 3L, integer(0)),
              class = c("bn", "parental")),
    structure(list(2L, 3L, integer(0)),
              class = c("bn", "parental")),
    structure(list(3L, 3L, integer(0)),
              class = c("bn", "parental")),
    structure(list(2:3, 3L, integer(0)),
              class = c("bn", "parental")),
    structure(list(integer(0), c(1L, 3L), integer(0)),
              class = c("bn", "parental")),
    structure(list(3L, c(1L, 3L), integer(0)),
              class = c("bn", "parental")),
    structure(list(integer(0), integer(0), 1L),
              class = c("bn", "parental")),
    structure(list(2L, integer(0), 1L),
              class = c("bn", "parental")),
    structure(list(integer(0), 1L, 1L),
              class = c("bn", "parental")),
    structure(list(integer(0), 3L, 1L),
              class = c("bn", "parental")),
    structure(list(integer(0), c(1L, 3L), 1L),
              class = c("bn", "parental")),
    structure(list(integer(0), integer(0), 2L),
              class = c("bn", "parental")),
    structure(list(2L, integer(0), 2L),
              class = c("bn", "parental")),
    structure(list(3L, integer(0), 2L),
              class = c("bn", "parental")),
    structure(list(2:3, integer(0), 2L),
              class = c("bn", "parental")),
    structure(list(integer(0), 1L, 2L),
              class = c("bn", "parental")),
    structure(list(integer(0), integer(0), 1:2),
              class = c("bn", "parental")),
    structure(list(2L, integer(0), 1:2),
              class = c("bn", "parental")),
    structure(list(integer(0), 1L, 1:2),
              class = c("bn", "parental"))),
    class = c("bn.list", "parental.list"))

  # each rows component refers to node1, node2, node3
  rows <- lapply(nets, optionsGivenGraph)

  groupScoresOld <- sapply(rows, getScoreFromRows)

  # sample group
  groupWeights <- exp(groupScoresOld - logsumexp(groupScoresOld))
  sampGroup <- sample.int(25, size = 1, prob = groupWeights)

  # sample 'node1' parents
  n1scoresGroup <- scoresParents[[node1]][rows[[sampGroup]][[1]]]
  n1probs <- exp(n1scoresGroup - logsumexp(n1scoresGroup))
  n1samp <- sample.int(length(n1scoresGroup), size = 1, prob = n1probs)

  # sample 'node2' parents
  n2scoresGroup <- scoresParents[[node2]][rows[[sampGroup]][[2]]]
  n2probs <- exp(n2scoresGroup - logsumexp(n2scoresGroup))
  n2samp <- sample.int(length(n2scoresGroup), size = 1, prob = n2probs)

  # sample 'node3' parents
  n3scoresGroup <- scoresParents[[node3]][rows[[sampGroup]][[3]]]
  n3probs <- exp(n3scoresGroup - logsumexp(n3scoresGroup))
  n3samp <- sample.int(length(n3scoresGroup), size = 1, prob = n3probs)

  # generate the new graph
  parents1 <- rows[[sampGroup]][[1]]
  new <- parentsTables[[node1]][parents1[n1samp], ]
  currentNetwork[[1]][[node1]] <- new[!is.na(new)]

  parents2 <- rows[[sampGroup]][[2]]
  new <- parentsTables[[node2]][parents2[n2samp], ]
  currentNetwork[[1]][[node2]] <- new[!is.na(new)]

  parents3 <- rows[[sampGroup]][[3]]
  new <- parentsTables[[node3]][parents3[n3samp], ]
  currentNetwork[[1]][[node3]] <- new[!is.na(new)]

  currentNetwork[[2]] <- routesAddEdges(currentNetwork[[2]],
                                        currentNetwork[[1]][[node1]],
                                        node1)
  currentNetwork[[2]] <- routesAddEdges(currentNetwork[[2]],
                                        currentNetwork[[1]][[node2]],
                                        node2)
  currentNetwork[[2]] <- routesAddEdges(currentNetwork[[2]],
                                        currentNetwork[[1]][[node3]],
                                        node3)
  currentNetwork[[4]][currentNetwork[[1]][[node1]], node1] <- 1
  if (length(currentNetwork[[1]][[node1]]) == 0){
    currentNetwork[[4]][, node1] <- 0
  } else {
    currentNetwork[[4]][-currentNetwork[[1]][[node1]], node1] <- 0
  }
  currentNetwork[[4]][currentNetwork[[1]][[node2]], node2] <- 1
  if (length(currentNetwork[[1]][[node2]]) == 0){
    currentNetwork[[4]][, node2] <- 0
  } else {
    currentNetwork[[4]][-currentNetwork[[1]][[node2]], node2] <- 0
  }
  currentNetwork[[4]][currentNetwork[[1]][[node3]], node3] <- 1
  if (length(currentNetwork[[1]][[node3]]) == 0){
    currentNetwork[[4]][, node3] <- 0
  } else {
    currentNetwork[[4]][-currentNetwork[[1]][[node3]], node3] <- 0
  }

  currentNetwork[[5]][node1] <- n1scoresGroup[n1samp]
  currentNetwork[[5]][node2] <- n2scoresGroup[n2samp]
  currentNetwork[[5]][node3] <- n3scoresGroup[n3samp]

  currentNetwork
}

#' Gibbs sampler for Bayesian Networks.
#'
#' Create a MCMC sampler for Bayesian Networks. The sampler samples Bayesian
#' Networks (ie models).
#'
#' @param data The data.
#' @param initial An object of class 'bn'. The starting value of the
#'                       MCMC.
#' @param prior A function that returns the prior score of the
#'                       supplied bn.
#' @param return Either "network" or "contingency".
#' @param logScoreFUN A list of four elements:
#'   \describe{
#'     \item{offline}{A function that computes the logScore of a Bayesian 
#'                    Network}
#'     \item{online}{A function that incrementally computes the logScore of a 
#'                   Bayesian Network}
#'     \item{local}{A function that computes the local logScore of a 
#'                  Bayesian Network}
#'     \item{prepare}{A function that prepares the data, and any further 
#'                    pre-computation required by the logScore functions.}
#'   }
#'   For Multinomial-Dirichlet models, \code{\link{logScoreMultDirFUN}} 
#'   returns the appropriate list; for Normal models with Zellner g-priors,
#'   \code{\link{logScoreNormalFUN}} returns the appropriate list.
#' @param logScoreParameters A list of parameters that are passed to
#'                       logScoreFUN.
#' @param constraint A matrix of dimension ncol(data) x ncol(data) giving
#'                       constraints to the sample space.
#'                       The (i, j) element is
#'                         1  if the edge i -> j is required
#'                         -1 if the edge i -> is excluded.
#'                         0  if the edge i -> j is not constrained.
#'                       The diagonal of constraint must be all 0.
#' @param statistics A named list of functions which should be applied to
#'   the current network after each step. Each function should accept an
#'   object of class \code{bn} and return a scalar output. Each item in
#'   the list must be named so that it can be referred to.
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node. A \code{NULL} value gives the default restriction 
#'   of 3.
#' @param moveprobs A numeric vector of length 3. Specifies the probability
#'   that moves updating the parent sets of 1, 2 and 3 nodes simultaneously.
#'   Must sum to 1.
#' @param verbose A logical of length 1, indicating whether verbose
#'                       output should be printed.
#' @param keepTape A logical of length 1, indicating whether a full log
#'                       ('tape') of the MCMC sampler should be kept.
#'                       Enabling this option can be very memory-intensive.
#' @param parentsTables A list of tables of the form returned by 
#'   \code{enumerateParentsTable()}
#' @param scoresParents A list of the form returned by 
#'   \code{scoreParentsTable()}
#' @return A function, which when called draws the next sample of the MCMC.
#' @export
#' @seealso \code{\link{BNSampler}}, \code{\link{BNSamplerBigFlips}}, 
#'   \code{\link{BNSamplerPT}}, \code{\link{BNSamplerMJ}},
#'   \code{\link{BNSamplerGrzeg}}. Internally uses 
#'   \code{\link{samplePair}} and \code{\link{sampleNode}}.
BNGibbsSampler <- function(data,
                           initial            = empty(ncol(data) - 1),
                           prior              = priorUniform(),
                           return             = "network",
                           logScoreFUN        = logScoreMultDirFUN(),
                           logScoreParameters = list(hyperparameters = "qi"),
                           constraint         = NULL,
                           statistics         = list(nEdges = nEdges),
                           maxNumberParents   = NULL,
                           moveprobs          = c(0.9, 0.1, 0),
                           verbose            = F,
                           keepTape           = F,
                           parentsTables      = NULL,
                           scoresParents      = NULL){
  stopifnot("bn" %in% class(initial),
            is.valid(initial),
            ncol(as.matrix(data)) ==   length(initial),
            is.function(prior),
            return               %in% c("network", "contingency"),
            class(statistics)     == "list",
            all(lapply(statistics, class) == "function"),
            all(nchar(names(statistics)) > 0),
            is.logical(keepTape),
            length(keepTape)      ==   1,
            sum(moveprobs)        ==   1)
  if (is.null(maxNumberParents)){
    maxNumberParents <- 3
  }

  numberOfNodes <- length(initial)
  nodesSeq <- seq_len(numberOfNodes)

  # Set up for fast computation of logScore
  logScoreLocalFUN <- logScoreFUN$local
  prepareDataFUN <- logScoreFUN$prepare
  logScoreParameters <- prepareDataFUN(data,
                                       logScoreParameters,
                                       checkInput = F)

  constraint <- setupConstraint(constraint, initial)
  required <- getRequiredFromConstraint(constraint)
  banned <- getBannedFromConstraint(constraint)

  if (is.null(parentsTables)){
    if (verbose){
      cat("Listing all possible parent sets\n")
      flush.console()
    }
    parentsTables <- enumerateParentsTable(numberOfNodes,
                                           maxNumberParents,
                                           required,
                                           banned,
                                           verbose = verbose)
  }
  if (is.null(scoresParents)){
    if (verbose){
      cat("Scoring all possible parent sets\n")
      flush.console()
    }
    scoresParents <- scoreParentsTable(parentsTables,
                                       logScoreLocalFUN,
                                       logScoreParameters,
                                       prior,
                                       verbose = verbose)
  }

  # The current MCMC state is stored a list of the form:
  # currentNetwork[[1]] is the bn
  # currentNetwork[[2]] is the routes matrix
  # currentNetwork[[3]] is the log prior
  currentNetwork <- vector(mode = "list", length = 5)
  currentNetwork[[1]] <- initial
  currentNetwork[[2]] <- routes(currentNetwork[[1]])
  currentNetwork[[3]] <- log(prior(currentNetwork[[1]]))
  currentNetwork[[4]] <- as.adjacency(currentNetwork[[1]])

  for (node in nodesSeq){
    localScore <- logScoreLocalFUN(node               = node,
                                   parents            = initial[[node]],
                                   logScoreParameters = logScoreParameters,
                                   cache              = new.env(),
                                   checkInput         = F)
    currentNetwork[[5]][node] <- localScore
  }

  if (!is.valid.prior(currentNetwork[[3]])){
    stop("Initial network has prior with 0 probability.")
  }

  rowsThatContain <- getRowsThatContain(numberOfNodes,
                                        parentsTables,
                                        maxNumberParents)

  allRows <- lapply(nodesSeq, function(node){
    seq_len(nrow(parentsTables[[node]]))
  })

  # Set up internal counters and logs etc
  nSteps <- 0

  if (return == "contingency"){
    count <- new.env(hash = T)
  }
  et <- matrix(0, nrow = numberOfNodes, ncol = numberOfNodes)
  etBinsIncrement <- 100
  etBinsSize <- 1000
  etbins <- matrix(ncol = numberOfNodes^2, nrow = etBinsIncrement)
  nBurnin <- 0

  statistics <- lapply(statistics, function(f){
    function(currentNetwork){
      f(currentNetwork[[1]])
    }
  })
  defaultStatistics <- list(logScores = function(currentNetwork){
    sum(currentNetwork[[5]])
  })
  statistics <- c(statistics, defaultStatistics)
  nStatistics <- length(statistics)
  statisticsTable <- matrix(ncol = nStatistics,
                            nrow = etBinsSize * etBinsIncrement)
  colnames(statisticsTable) <- names(statistics)

  if (isTRUE(keepTape)){
    tapeSizeIncrement <- 500000
    tapeColumns <- c("movetype", "logAccProb", "accepted")
    numberTapeColumns <- length(tapeColumns)
    tape <- matrix(nrow = 0, ncol = numberTapeColumns)
    tapeProposals <- character(length = 0)
    colnames(tape) <- tapeColumns
  }

  returnDiagnostics <- function(){
    list(nSteps = nSteps)
  }

  returnTape <- function(){
    if (isTRUE(keepTape)){
      data.frame(tape[seq_len(nSteps), ],
                 proposals = tapeProposals[seq_len(nSteps)])
    }
    else {
      warning("Tape not kept for this MCMC run")
      data.frame()
    }
  }

  lengthenTape <- function(){
    if (nSteps %% tapeSizeIncrement == 0){
      tapeTemp <- tape
      tape <<- matrix(nrow = nSteps + tapeSizeIncrement,
                      ncol = numberTapeColumns)
      colnames(tape) <<- tapeColumns
      tape[seq_len(nSteps), ] <<- tapeTemp

      tapeProposalsTemp <- tapeProposals
      tapeProposals <<- character(length = nSteps + tapeSizeIncrement)
      tapeProposals[seq_len(nSteps)] <<- tapeProposalsTemp
    }
  }

  updateTape <- function(nSteps, currentNetwork){
    tape[nSteps, 1] <<- -1
    tape[nSteps, 2] <<- -1
    tape[nSteps, 3] <<- 1
    tapeProposals[nSteps] <<- as.character(currentNetwork[[1]], pretty = T)
  }

  updateET <- function(currentNetwork, nSteps, nBurnin){
    if (nSteps > nBurnin){
      et <<- et + currentNetwork[[4]]
      lengthenETBins(nSteps, nBurnin)
      row <- ((nSteps - nBurnin - 1) %/% etBinsSize) + 1
      etbins[row, ] <<- as.vector(et)
      if ((nSteps - nBurnin) %% etBinsSize == 0){
        et <<- matrix(0, nrow = numberOfNodes, ncol = numberOfNodes)
      }
    }
  }

  lengthenETBins <- function(nSteps, nBurnin){
    if ((nSteps - nBurnin) %% (etBinsSize * etBinsIncrement) == 0){
      temp <- etbins
      nRowsPrev <- nrow(etbins)
      etbins <<- matrix(nrow = nRowsPrev + etBinsIncrement,
                        ncol = numberOfNodes^2)
      etbins[seq_len(nRowsPrev), ] <<- temp
    }
  }

  updateStatistics <- function(currentNetwork, nSteps, nBurnin){
    lengthenStatistics(nSteps, nBurnin)
    if (nSteps > nBurnin){
      step <- nSteps - nBurnin
      for (i in seq_along(statistics)){
        statisticsTable[step, i] <<- statistics[[i]](currentNetwork)
      }
    }
  }

  lengthenStatistics <- function(nSteps, nBurnin){
    if ((nSteps - nBurnin) %% (etBinsSize * etBinsIncrement) == 0){
      temp <- statisticsTable
      nRowsPrev <- nrow(temp)
      nRowsNew <- nRowsPrev + etBinsSize * etBinsIncrement
      statisticsTable <<- matrix(ncol = nStatistics, nrow = nRowsNew)
      colnames(statisticsTable) <<- names(statistics)
      statisticsTable[seq_len(nRowsPrev), ] <<- temp
    }
  }

  sampler <- function(x,
                      verbose = F,
                      returnDiagnostics = F,
                      debugAcceptance = F,
                      returnTape = F,
                      burnin = 0){
    if (isTRUE(returnDiagnostics)) return(returnDiagnostics())
    if (isTRUE(returnTape)) return(returnTape())
    if (isTRUE(debugAcceptance)) browser()
    if (isTRUE(keepTape)) lengthenTape()

    nBurnin <<- burnin
    nSteps <<- nSteps + 1

    u <- runif(1, min = 0, max = 1)
    if (u < moveprobs[1]){
      currentNetwork <<- sampleNode(currentNetwork  = currentNetwork,
                                    numberOfNodes   = numberOfNodes,
                                    nodesSeq        = nodesSeq,
                                    scoresParents   = scoresParents,
                                    parentsTables   = parentsTables,
                                    allRows         = allRows,
                                    rowsThatContain = rowsThatContain)
    } else if (u < moveprobs[1] + moveprobs[2]){
      currentNetwork <<- samplePair(currentNetwork  = currentNetwork,
                                    numberOfNodes   = numberOfNodes,
                                    nodesSeq        = nodesSeq,
                                    scoresParents   = scoresParents,
                                    parentsTables   = parentsTables,
                                    allRows         = allRows,
                                    rowsThatContain = rowsThatContain)
    } else if (u < sum(moveprobs[1:3])) {
      currentNetwork <<- sampleTriple(currentNetwork  = currentNetwork,
                                      numberOfNodes   = numberOfNodes,
                                      nodesSeq        = nodesSeq,
                                      scoresParents   = scoresParents,
                                      parentsTables   = parentsTables,
                                      allRows         = allRows,
                                      rowsThatContain = rowsThatContain,
                                      logScoreFUN     = logScoreFUN,
                                      logScoreParameters = logScoreParameters)
    } else {
      stop("Not implemented")
    }

    if (isTRUE(keepTape)) updateTape(nSteps, currentNetwork)
    if (isTRUE(debugAcceptance)) browser()

    updateET(currentNetwork, nSteps, nBurnin)
    updateStatistics(currentNetwork, nSteps, nBurnin)

    if (return == "network"){
      currentNetwork[[1]]
    } else if (return == "contingency") {
      if (nSteps > nBurnin){
        id <- as.character(currentNetwork[[1]])
        if (is.null(count[[id]])){
          count[[id]] <<- 1L
        }
        else {
          count[[id]] <<- count[[id]] + 1L
        }
        currentNetwork[[1]]
      } else {
        NULL
      }
    }
  }
  class(sampler) <- c("sampler", "function")
  sampler
}
