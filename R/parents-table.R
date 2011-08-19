# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Table of all possible parents of a node.
#'
#' Creates a list of tables with rows corresponding to the possible parent 
#' sets of each node.
#'
#' @param numberOfNodes Number of nodes
#' @param maxNumberParents The maximum indegree of the node.
#' @param required A list of numeric vectors. Each component is a numeric 
#'   vector containing the nodes that must be included in every parent set 
#'   for this node.
#' @param banned A list of numeric vectors. Each component is a numeric 
#'   vector containing the nodes that must be excluded from every parent 
#'   set for this node.
#' @param verbose A logical of length 1, indicating whether verbose output
#'   should be printed.
#' @return A list of matrices of the form returned by 
#'   \code{enumerateParentsTableNode}.
#' @export
#' @seealso \code{\link{enumerateParentsTableNode}},
#'   \code{\link{scoreParentsTable}}, \code{\link{whichParentSetRows}},
#'   \code{\link{getRowsThatContain}}
enumerateParentsTable <- function(numberOfNodes,
                                  maxNumberParents,
                                  required,
                                  banned,
                                  verbose = F){
  nodesSeq <- seq_len(numberOfNodes)
  
  if (isTRUE(verbose)){
    progress <- txtProgressBar(max = numberOfNodes, style = 3)
    setTxtProgressBar(progress, 0)
  }
  
  parentsTables <- vector("list", numberOfNodes)
  for (node in nodesSeq){
    parentsTables[[node]] <- enumerateParentsTableNode(
                               node             = node,
                               numberOfNodes    = numberOfNodes,
                               maxNumberParents = maxNumberParents,
                               required         = required[[node]],
                               banned           = banned[[node]])
    if (isTRUE(verbose)){
      setTxtProgressBar(progress, node)
    }
  }
  if (isTRUE(verbose)){
    close(progress)
  }
  parentsTables
}

#' Table of all possible parents of a node.
#'
#' Creates a matrix, with each row being a parent set.
#' The is is created subject to the supplied indegree restriction, and the 
#' the supplied \code{required} and \code{banned} restrictions.
#'
#' @param node A node
#' @param numberOfNodes Number of nodes
#' @param maxNumberParents The maximum indegree of the node.
#' @param required A numeric vector. The nodes that must be included in 
#'   every parent set for this node.
#' @param banned A numeric vector. The nodes that must be excluded from 
#'   every parent set for this node.
#' @return A matrix with \code{maxNumberParents} columns.
#'   Each row is a possible parent set for node \code{node}, accounting for 
#'   the restrictions given by \code{required} and \code{banned}.
#'   Entries that are \code{NA} indicate no parent. e.g. there is only one 
#'   parent, the other entries will be \code{NA}.
#' @export
#' @seealso \code{\link{enumerateParentsTable}}
enumerateParentsTableNode <- function(node,
                                      numberOfNodes,
                                      maxNumberParents,
                                      required,
                                      banned){
  nodesSeq <- seq_len(numberOfNodes)
  nonDescendants <- setdiff3(nodesSeq, node)
  nonDescendants <- setdiff3(nonDescendants, banned)
  
  padWithNAs <- function(x){
    lenx <- length(x)
    if (lenx < maxNumberParents){
      c(x, rep(NA, maxNumberParents - lenx))
    } else {
      x
    }
  }
  
  parents <- enumerateParents(nonDescendants,
                              maxNumberParents = maxNumberParents,
                              required        = required)
  parents <- lapply(parents, padWithNAs)
  
  out <- matrix(unlist(parents), ncol = maxNumberParents, byrow = T)
  storage.mode(out) <- "integer"
  out
}

#' Score a node-level parentsTable.
#'
#' Computes the scores of all the Bayesian Networks, with parent sets 
#' corresponding to each row of a single component of a \code{parentsTable}.
#'
#' @param parentsTables A component of a \code{parentsTables}, of the form 
#'   created by \code{enumerateParentsTable()}.
#' @param logScoreLocalFUN A function that computes the local logScore 
#'   of a Bayesian Network.
#' @param logScoreParameters A list of parameters that are passed to
#'   \code{logScoreFUN}.
#' @param loalPriors A list of functions of the same length as the number of
#'   nodes that returns the local prior score of the corresponding node, given a
#'   numeric vector of parents.
#' @param verbose A logical of length 1, indicating whether verbose output
#'   should be printed.
#' @return List of numeric vectors of scores.
#' @export
#' @seealso \code{\link{scoreParentsTableNode}},
#'   \code{\link{enumerateParentsTable}}, \code{\link{whichParentSetRows}},
#'   \code{\link{getRowsThatContain}}
scoreParentsTable <- function(parentsTables,
                              logScoreLocalFUN,
                              logScoreParameters,
                              prior,
                              verbose = F){
  numberOfNodes <- length(parentsTables)
  numberOfRows <- sapply(parentsTables, nrow)
  totalNumberRows <- sum(numberOfRows)
  nodesSeq <- seq_len(numberOfNodes)
  if (isTRUE(verbose)){
    progress <- txtProgressBar(max = totalNumberRows, style = 3)
    setTxtProgressBar(progress, 0)
  }

  scoresList <- vector("list", numberOfNodes)
  for (node in nodesSeq){
    localPrior <- prior[[node]]
    scoresList[[node]] <- scoreParentsTableNode(
                          node               = node,
                          parentsTable       = parentsTables[[node]],
                          logScoreLocalFUN   = logScoreLocalFUN,
                          logScoreParameters = logScoreParameters,
                          localPrior         = localPrior)
    if (isTRUE(verbose)){
      new <- getTxtProgressBar(progress) + numberOfRows[[node]]
      setTxtProgressBar(progress, new)
    }
  }
  if (isTRUE(verbose)){
    close(progress)
  }
  scoresList
}

#' Score a node-level parentsTable.
#'
#' Computes the local scores of all the Bayesian Networks, with node 
#' \code{node} set to parent sets corresponding to each row of a single 
#' component of a \code{parentsTable}.
#'
#' @param node A node. A numeric vector of length 1.
#' @param parentsTable A component of a \code{parentsTable}, of the form 
#'   created by \code{enumerateParentsTable()}.
#' @param logScoreLocalFUN A function that computes the local logScore
#'   of a Bayesian Network.
#' @param logScoreParameters A list of parameters that are passed to
#'   \code{logScoreFUN}.
#' @param localPrior  function that returns the prior score of the
#'   supplied parents.
#' @return A numeric vector of scores.
#' @export
#' @seealso \code{\link{scoreParentsTable}}
scoreParentsTableNode <- function(node,
                                  parentsTable,
                                  logScoreLocalFUN,
                                  logScoreParameters,
                                  localPrior){
  i <- 0
  nr <- nrow(parentsTable)
  scores <- vector("numeric", length = nr)
  while (i <= nr){
    parents <- parentsTable[i, ]
    parents <- parents[!is.na(parents)]
    scores[i] <- logScoreLocalFUN(node               = node,
                                  parents            = parents,
                                  logScoreParameters = logScoreParameters,
                                  cache              = new.env(),
                                  checkInput         = F)
    scores[i] <- log(localPrior(parents)) + scores[i]
    i <- i + 1
  }
  scores
}

#' Create lookup table for parentsTable.
#'
#' Creates a list that allows quick lookup of a parentsTable. This is 
#' needed for \code{whichParentSetRows}.
#'
#' @param numberOfNodes The number of nodes. A numeric vector.
#' @param parentsTables A list of tables of the form returned by 
#'   \code{enumerateParentsTable()}
#' @param maxNumberParents The maximum indegree for each node. A numeric 
#'   vector of length 1.
#' @return A list of length \code{numberOfNodes}. Each component of this 
#'   list is a list of length \code{numberOfNodes}, the i^{th} component of 
#'   which is a numeric vector containing the indicies of 
#'   \code{parentsTables[[node]]} that correspond to parent sets including 
#'   \code{i}.
#' @export
#' @seealso \code{\link{enumerateParentsTable}},
#'   \code{\link{scoreParentsTable}}
getRowsThatContain <- function(numberOfNodes,
                               parentsTables,
                               maxNumberParents){
  nodesSeq <- seq_len(numberOfNodes)
  lapply(nodesSeq, function(node){
    lapply(nodesSeq, function(parent){
      size <- nrow(parentsTables[[node]])
      ind <- match(as.vector(parentsTables[[node]]), parent, 0, NULL)
      tab <- matrix2(ind, nrow = size, ncol = maxNumberParents)
      tab <- rowSums2(tab)
      .Internal(which(tab > 0))
    })
  })
}

#' Find relevants rows of a parentsTable.
#'
#' Finds the rows of a parentsTable that correspond to parent sets that 
#' could be added as parents of node \code{node}, given some set of 
#' nodes \code{nonDescendants} that can be added as parents without 
#' creating a cycle in the graph.
#'
#' Note that nodes that are banned do not need to be accounted 
#' for in the \code{nonDescendants} argument, since these should be 
#' accounted for when the parentsTable is created. Required nodes must be 
#' included in \code{nonDescendants}.
#'
#' @param node The node. A numeric vector of length 1.
#' @param nonDescendants The nodes that can be added as descendants of 
#'   \code{node}. A numeric vector.
#' @param needOneOf Nodes the MUST be included as parents of the node. This
#'   is NOT the same as the required/banned list; this is to allow the
#'   space to be partitioned.
#' @param numberOfNodes The number of nodes in the network. A numeric vector 
#'   of length 1.
#' @param allRows The vector 1:nrow(parentsTables). (Supplied as an 
#'   argument for possible speed gain)
#' @param rowsThatContain A list of the form created by 
#'   \code{getRowsThatContain()}
#' @return A numeric vector.
#' @export
#' @seealso \code{\link{enumerateParentsTable}},
#'   \code{\link{scoreParentsTable}}
whichParentSetRows <- function(node,
                               nonDescendants,
                               needOneOf = NULL,
                               numberOfNodes,
                               allRows,
                               rowsThatContain){
  nodesSeq <- seq_len(numberOfNodes)
  nodesNotAllowed <- setdiff3(nodesSeq, nonDescendants)
  rowsNotAllowed <- rowsThatContain[[node]][nodesNotAllowed]
  rowsNotAllowed <- unlist(rowsNotAllowed, use.names = F)
  rowsNeeded <- allRows[[node]]
  if (!is.null(needOneOf) && length(needOneOf) > 0){
    if (is.list(needOneOf) && length(needOneOf) == 1){
      needOneOf <- needOneOf[[1]]
    }
    if (is.list(needOneOf)){
      rowsNeeded <- lapply(needOneOf, function(needed){
        unlist(rowsThatContain[[node]][needed], use.names = F)
      })
      if (length(rowsNeeded) == 2){
        rowsNeeded <- intersect2(rowsNeeded[[1]], rowsNeeded[[2]])
      } else {
        rowsNeeded <- do.call("intersection", rowsNeeded)
      }
    } else {
      rowsNeeded <- rowsThatContain[[node]][needOneOf]
      rowsNeeded <- unlist(rowsNeeded, use.names = F)
      rowsNeeded <- .Internal(unique(rowsNeeded, F, F))
    }
  }
  setdiff3(rowsNeeded, rowsNotAllowed)
}
