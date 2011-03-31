# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Sample the parents of a single node
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
  currentNetwork
}

#' Sample the parents of a pair of nodes
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

  rows1 <- vector("list", 2)
  rows2 <- vector("list", 2)

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
  newNonDescendants2 <- nonDescendants2

  rows2[[1]] <- whichParentSetRows(node            = node2,
                                   nonDescendants  = newNonDescendants2,
                                   numberOfNodes   = numberOfNodes,
                                   allRows         = allRows,
                                   rowsThatContain = rowsThatContain)
  group1Score <- logsumexp(scoresParents[[node1]][rows1[[1]]]) +
                 logsumexp(scoresParents[[node2]][rows2[[1]]])

  # group2
  newNonDescendants1 <- intersect2(nonDescendants1, descendants2)

  # get the conditional probability for the parents of
  # node 'node1', given the rest of the graph
  rows1[[2]] <- whichParentSetRows(node            = node1,
                                   nonDescendants  = newNonDescendants1,
                                   numberOfNodes   = numberOfNodes,
                                   allRows         = allRows,
                                   rowsThatContain = rowsThatContain)

  # haveNewDescendants == T
  # nonDescendants2 = nonDescendants2 + descendants1

  newNonDescendants2 <- setdiff3(nonDescendants2, descendants1)
  rows2[[2]] <- whichParentSetRows(node            = node2,
                                   nonDescendants  = newNonDescendants2,
                                   numberOfNodes   = numberOfNodes,
                                   allRows         = allRows,
                                   rowsThatContain = rowsThatContain)

  group2Score <- logsumexp(scoresParents[[node1]][rows1[[2]]]) +
                 logsumexp(scoresParents[[node2]][rows2[[2]]])

  # sample group
  groupScoresOld <- c(group1Score, group2Score)
  groupWeights <- exp(groupScoresOld - logsumexp(groupScoresOld))
  n2SampGroup <- sample.int(2, size = 1, prob = groupWeights)
  
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
  currentNetwork
}

#' Gibbs sampler for Bayesian Networks
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
#' @param logScoreFUN A list of three elements:
#'                         offline: A function that computes the logScore
#'                                  of a Bayesian Network
#'                         online:  A function that incrementally computes
#'                                  the logScore of a Bayesian Network
#'                         prepare: A function that prepares the data, and
#'                                  any further pre-computation required by
#'                                  the logScore functions.
#' @param logScoreParameters A list of parameters that are passed to
#'                       logScoreFUN.
#' @param constraint A matrix of dimension ncol(data) x ncol(data) giving
#'                       constraints to the sample space.
#'                       The (i, j) element is
#'                         1  if the edge i -> j is required
#'                         -1 if the edge i -> is excluded.
#'                         0  if the edge i -> j is not constrained.
#'                       The diagonal of constraint must be all 0.
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node.
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
BNGibbsSampler <- function(data,
                           initial,
                           prior,
                           return      = "network",
                           logScoreFUN = defaultLogScoreFUN(),
                           logScoreParameters = list(hyperparameters = "qi"),
                           constraint  = NULL,
                           maxNumberParents,
                           moveprobs = c(0.9, 0.1, 0),
                           verbose     = F,
                           keepTape    = F,
                           parentsTables,
                           scoresParents){
  stopifnot("bn" %in% class(initial),
            is.valid(initial),
            ncol(as.matrix(data)) ==   length(initial),
            is.function(prior),
            return               %in% c("network", "contingency", "sha256"),
            is.logical(keepTape),
            length(keepTape)      ==   1,
            sum(moveprobs)        ==   1)

  numberOfNodes <- length(initial)
  nodesSeq <- seq_len(numberOfNodes)

  # Set up for fast computation of logScore
  logScoreOfflineFUN <- logScoreFUN$offline
  prepareDataFUN <- logScoreFUN$prepare
  logScoreParameters <- prepareDataFUN(data,
                                       logScoreParameters,
                                       checkInput = F)

  constraint <- setupConstraint(constraint, initial)
  required <- getRequiredFromConstraint(constraint)
  banned <- getBannedFromConstraint(constraint)
  
  if (missing(parentsTables)){
    parentsTables <- enumerateParentsTable(numberOfNodes,
                                           maxNumberParents,
                                           required,
                                           banned,
                                           verbose = verbose)
  }
  if (missing(scoresParents)){
    scoresParents <- scoreParentsTable(parentsTables,
                                       logScoreOfflineFUN,
                                       logScoreParameters,
                                       prior,
                                       verbose = verbose)
  }
  

  # The current MCMC state is stored a list of the form:
  # currentNetwork[[1]] is the bn
  # currentNetwork[[2]] is the routes matrix
  # currentNetwork[[3]] is the log prior
  currentNetwork <- vector(mode = "list", length = 3)
  currentNetwork[[1]] <- initial
  currentNetwork[[2]] <- routes(currentNetwork[[1]])
  currentNetwork[[3]] <- log(prior(currentNetwork[[1]]))
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
    lookup <- new.env(hash = T)
  } else if (return == "sha256"){
    lookup <- new.env(hash = T)
  }

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

  function(x,
           verbose = F,
           returnDiagnostics = F,
           debugAcceptance = F,
           returnTape = F,
           burnin = 0){
    if (isTRUE(returnDiagnostics)) return(returnDiagnostics())
    if (isTRUE(returnTape)) return(returnTape())
    if (isTRUE(debugAcceptance)) browser()
    if (isTRUE(keepTape)) lengthenTape()

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
    } else {
      stop("Not implemented")
    }

    if (isTRUE(keepTape)) updateTape(nSteps, currentNetwork)
    if (isTRUE(debugAcceptance)) browser()

    if (return == "network"){
      currentNetwork[[1]]
    } else if (return == "contingency") {
      if (nSteps > burnin){
        id <- digest(currentNetwork[[1]], algo = "sha256")
        if (is.null(count[[id]])){
          count[[id]] <<- 1L
          lookup[[id]] <<- currentNetwork[[1]]
        }
        else {
          count[[id]] <<- count[[id]] + 1L
        }
        currentNetwork[[1]]
      } else {
        NULL
      }
    } else if (return == "sha256"){
      id <- digest(currentNetwork[[1]], algo = "sha256")
      if (is.null(count[[id]])){
        lookup[[id]] <<- currentNetwork[[1]]
      } else {
        lookup[[id]] <<- currentNetwork[[1]]
      }
      id
    }
  }
}
