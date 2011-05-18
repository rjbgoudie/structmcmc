# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' (Log) Number of neighbouring networks.
#'
#' Returns the number of acyclic graphs that can be formed by adding, 
#' removing or flipping a single edge of the current network
#'
#' @param routes The routes matrix of the network
#' @param adjacency The adjacency matrix of the network
#' @param constraintT The transpose of a constraint matrix
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node.
#' @return A numeric of length 1. The log number of neighbouring graphs.
#' @export
#' @seealso \code{\link{BNSampler}}
logNumMHNeighbours <- function(routes,
                               adjacency,
                               constraintT,
                               maxNumberParents = Inf){
  removable <- transposeEdgeIsRemovable(routes      = routes,
                                        adjacency   = adjacency,
                                        constraintT = constraintT)
  addable <- transposeEdgeIsAddable(routes           = routes,
                                    adjacency        = adjacency,
                                    constraintT      = constraintT,
                                    maxNumberParents = maxNumberParents)
  flippable <- edgeIsFlippable(routes           = routes,
                               adjacency        = adjacency,
                               constraintT      = constraintT,
                               maxNumberParents = maxNumberParents)
  # Note that flippable is the edge itself, not the transpose
  # So an edge i -> j that is flippable will also count under j -> i
  # for removable.
  log(length(adjacency[flippable | addable | removable]))
}

#' Find togglable edges.
#'
#' Finds "edge-locations" in the graph that can be added or removed without
#' introducing a cycle into the graph
#'
#' @param routes The routes matrix of the network
#' @param adjacency The adjacency matrix of the network
#' @param constraintT The transpose of a constraint matrix
#' @return A logical matrix of the same dimension as the supplied matrices, 
#'   with entries indicating whether the corresponding edge can be added or 
#'   removed without introducing a cycle. NOTE THIS IS TRANSPOSE OF EXPECTED
#' @export
#' @seealso \code{\link{BNSampler}}, \code{\link{transposeEdgeIsAddable}},
#'   \code{\link{transposeEdgeIsTogglable}}, \code{\link{edgeIsFlippable}}
transposeEdgeIsRemovable <- function(routes, adjacency, constraintT){
  routes == 0 & t(adjacency) == 1 & constraintT == 0
}

#' Find togglable edges.
#'
#' Finds "edge-locations" in the graph that can be added or removed without
#' introducing a cycle into the graph
#'
#' @param routes The routes matrix of the network
#' @param adjacency The adjacency matrix of the network
#' @param constraintT The transpose of a constraint matrix
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node.
#' @return A logical matrix of the same dimension as the supplied matrices, 
#'   with entries indicating whether the corresponding edge can be added or 
#'   removed without introducing a cycle. NOTE THIS IS TRANSPOSE OF EXPECTED
#' @export
#' @seealso \code{\link{BNSampler}}, \code{\link{transposeEdgeIsRemovable}},
#'   \code{\link{transposeEdgeIsTogglable}}, \code{\link{edgeIsFlippable}}
transposeEdgeIsAddable <- function(routes,
                                   adjacency,
                                   constraintT,
                                   maxNumberParents){
  addable <- routes == 0 & t(adjacency) == 0 & constraintT == 0
  canAddMoreParents <- colSums(adjacency) < maxNumberParents
  addable[!canAddMoreParents, ] <- F
  addable
}

#' Find togglable edges.
#'
#' Finds "edge-locations" in the graph that can be added or removed without
#' introducing a cycle into the graph
#'
#' @param routes The routes matrix of the network
#' @param adjacency The adjacency matrix of the network
#' @param constraintT The transpose of a constraint matrix
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node.
#' @return A logical matrix of the same dimension as the supplied matrices, 
#'   with entries indicating whether the corresponding edge can be added or 
#'   removed without introducing a cycle. NOTE THIS IS TRANSPOSE OF EXPECTED
#' @export
#' @seealso \code{\link{BNSampler}}, \code{\link{transposeEdgeIsAddable}},
#'   \code{\link{transposeEdgeIsRemovable}}, \code{\link{edgeIsFlippable}}
transposeEdgeIsTogglable <- function(routes,
                                     adjacency,
                                     constraintT,
                                     maxNumberParents = Inf){
  removable <- transposeEdgeIsRemovable(routes      = routes,
                                        adjacency   = adjacency,
                                        constraintT = constraintT)
  addable <- transposeEdgeIsAddable(routes           = routes,
                                    adjacency        = adjacency,
                                    constraintT      = constraintT,
                                    maxNumberParents = maxNumberParents)
  removable + addable == T
}

#' Find flippable edges.
#'
#' Finds edges in the graph whose direction can be reversed ("flipped")
#' without introducing a cycle into the graph.
#'
#' @param routes The routes matrix of the network
#' @param adjacency The adjacency matrix of the network
#' @param constraintT The transpose of a constraint matrix
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node.
#' @return A logical matrix of the same dimension as the supplied matrices, 
#'   with entries indicating whether the corresponding edge can be flipped
#'   without introducing a cycle. NOTE THIS IS TRANSPOSE OF EXPECTED
#' @export
#' @seealso \code{\link{BNSampler}}, \code{\link{transposeEdgeIsAddable}},
#'   \code{\link{transposeEdgeIsRemovable}},
#'   \code{\link{transposeEdgeIsTogglable}}
edgeIsFlippable <- function(routes, adjacency, constraintT, maxNumberParents){
  flippable <- routes == 1 & adjacency == 1 & constraintT == 0
  canAddMoreParents <- colSums(adjacency) < maxNumberParents
  flippable[!canAddMoreParents, ] <- F
  flippable
}

#' Create a MCMC sampler (MC^3) for Bayesian Networks.
#' 
#' The sampler samples Bayesian Networks (ie models).
#'
#' @param data The data.
#' @param initial An object of class 'bn'. The starting value of the MCMC.
#' @param prior A function that returns the prior score of the supplied bn.
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
#'   \code{\link{logScoreZellnerFUN}} returns the appropriate list.
#' @param logScoreParameters A list of parameters that are passed to
#'   logScoreFUN.
#' @param constraint A matrix of dimension ncol(data) x ncol(data) giving
#'   constraints to the sample space. The (i, j) element is:
#'   \describe{
#'     \item{1}{if the edge i -> j is required}
#'     \item{-1}{if the edge i -> is excluded.}
#'     \item{0}{if the edge i -> j is not constrained.}
#'   }
#'   The diagonal of constraint must be all 0.
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node. The default value, which is used for \code{NULL}
#'   is to not constrain the maximum indegree, i.e. to use
#'   \code{ncol(data) - 1}.
#' @param verbose A logical of length 1, indicating whether verbose
#'  output should be printed.
#' @param keepTape A logical of length 1, indicating whether a full log
#'   (\code{tape}) of the MCMC sampler should be kept. Enabling this option 
#'   can be very memory-intensive.
#' @return A function, which when called draws the next sample of the MCMC.
#' @export
#' @seealso \code{\link{draw}}. \code{\link{BNSamplerMJ}},
#'   \code{\link{BNSamplerPT}}, \code{\link{BNGibbsSampler}}. Example priors
#'   \code{\link{priorGraph}}, \code{\link{priorUniform}}
#' @examples
#' x1 <- factor(c("a", "a", "g", "c", "c", "a", "g", "a", "a"))
#' x2 <- factor(c(2, 2, 4, 3, 1, 4, 4, 4, 1))
#' x3 <- factor(c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
#' x <- data.frame(x1 = x1, x2 = x2, x3 = x3)
#' 
#' initial <- empty(3, "bn")
#' prior <- priorUniform()
#' 
#' sampler <- BNSampler(data = x, initial = initial, prior = prior)
#' samples <- draw(sampler, n = 100, burnin = 10)
#' 
#' x <- bnpostmcmc(sampler, samples)
#' ep(x)
BNSampler <- function(data,
                      initial,
                      prior,
                      return      = "network",
                      logScoreFUN = logScoreMultDirFUN(),
                      logScoreParameters = list(hyperparameters = "qi"),
                      constraint  = NULL,
                      maxNumberParents = NULL,
                      verbose     = F,
                      keepTape    = F){
  if (is.null(maxNumberParents)){
    maxNumberParents <- ncol(data) - 1
  }
  stopifnot("bn"                  %in% class(initial),
            is.valid(initial),
            ncol(as.matrix(data)) ==   length(initial),
            is.function(prior),
            return                %in% c("network", "contingency", "sha256"),
            inherits(maxNumberParents, "numeric") ||
              inherits(maxNumberParents, "integer"),
            is.logical(keepTape),
            length(keepTape)      ==   1)

  numberOfNodes <- length(initial)
  nodesSeq <- seq_len(numberOfNodes)
  cache <- new.env(hash = T, size = 10000L)
  
  # Set up for fast computation of logScore
  logScoreOfflineFUN <- logScoreFUN$offline
  logScoreOnlineFUN <- logScoreFUN$online
  prepareDataFUN <- logScoreFUN$prepare
  logScoreParameters <- prepareDataFUN(data,
                                       logScoreParameters,
                                       checkInput = F)

  constraint <- setupConstraint(constraint, initial)
  constraintT <- t(constraint)

  # The current MCMC state is stored a list of the form:
  # currentNetwork[[1]] is the bn
  # currentNetwork[[2]] is the routes matrix
  # currentNetwork[[3]] is the log prior
  # currentNetwork[[4]] is the adjacency matrix
  # currentNetwork[[5]] is the log number of neighbours
  currentNetwork <- list(5, mode = "list")
  currentNetwork[[1]] <- initial
  currentScore <- logScoreOfflineFUN(x                  = currentNetwork[[1]],
                                     logScoreParameters = logScoreParameters,
                                     cache              = cache)
  currentNetwork[[2]] <- routes(currentNetwork[[1]])
  currentNetwork[[3]] <- log(prior(currentNetwork[[1]]))
  if (!is.valid.prior(currentNetwork[[3]])){
    stop("Initial network has prior with 0 probability.")
  }
  currentNetwork[[4]] <- as.adjacency(currentNetwork[[1]])
  currentNetwork[[5]] <- logNumMHNeighbours(routes      = currentNetwork[[2]],
                                            adjacency   = currentNetwork[[4]],
                                            constraintT = constraintT)

  # Set up internal counters and logs etc
  nAccepted <- 0
  nSteps <- 0
  nMHProposals <- 0
  nMHAccepted <- 0
  nCurrentGraphIsAMode <- 0
  if (return == "contingency"){
    count <- new.env(hash = T)
    lookup <- new.env(hash = T)
  } else if (return == "sha256"){
    lookup <- new.env(hash = T)
  }

  if (isTRUE(keepTape)){
    tapeSizeIncrement <- 500000
    tapeColumns <- c("movetype", "logAccProb", "accepted",
                     "currIsMode", "propIsMode")
    numberTapeColumns <- length(tapeColumns)
    tape <- matrix(nrow = 0, ncol = numberTapeColumns)
    tapeProposals <- character(length = 0)
    colnames(tape) <- tapeColumns
  }

  returnDiagnostics <- function(){
    list(nAccepted            = nAccepted,
         nSteps               = nSteps,
         acceptanceRate       = nAccepted/nSteps,
         nMHProposals         = nMHProposals,
         nMHAccepted          = nMHAccepted,
         MHAcceptanceRate     = nMHAccepted/nMHProposals,
         nCurrentGraphIsAMode = nCurrentGraphIsAMode)
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

  updateTape <- function(nSteps,
                         currentNetwork,
                         movetype,
                         logAccProb,
                         accepted){
    tape[nSteps, 1] <<- if (movetype == "mh") 0 else 1
    tape[nSteps, 2] <<- logAccProb
    tape[nSteps, 3] <<- accepted
    tapeProposals[nSteps] <<- as.character(currentNetwork[[1]], pretty = T)
  }

  function(x,
           verbose           = F,
           returnDiagnostics = F,
           debugAcceptance   = F,
           returnTape        = F,
           burnin            = 0){
    if (isTRUE(returnDiagnostics)) return(returnDiagnostics())
    if (isTRUE(returnTape)) return(returnTape())
    if (isTRUE(debugAcceptance)) browser()
    if (isTRUE(keepTape)) lengthenTape()

    nSteps <<- nSteps + 1

    proposalNetwork <- currentNetwork

    # function for generating a proposal
    # returns the acceptance probability

    logp <- log(runif(1, min = 0, max = 1))

    # propose mh move
    nMHProposals <<- nMHProposals + 1

    # count the number of proposals and select one
    canAddOrRemove <- transposeEdgeIsTogglable(
                        routes           = currentNetwork[[2]],
                        adjacency        = currentNetwork[[4]],
                        constraintT      = constraintT,
                        maxNumberParents = maxNumberParents)
    canFlip <- edgeIsFlippable(routes           = currentNetwork[[2]],
                               adjacency        = currentNetwork[[4]],
                               constraintT      = constraintT,
                               maxNumberParents = maxNumberParents)

    nonCycleInducing <- which(canAddOrRemove, arr.ind = T)
    nonCycleInducingFlips <- which(canFlip, arr.ind = T)

    nNonCycleInducing <- nrow(nonCycleInducing)
    nNonCycleInducingFlips <- nrow(nonCycleInducingFlips)
    numberOfMoves <- nNonCycleInducing + nNonCycleInducingFlips
    if (numberOfMoves < 1){
      stop("No available moves")
    }
    select <- sample.int(numberOfMoves, size = 1)

    if (select <= nNonCycleInducing){
      # swapped to save transposing the matrix
      j <- nonCycleInducing[[select, 1]]
      i <- nonCycleInducing[[select, 2]]

      # the condition is a marginally faster i %in% currentNetwork[[1]][[j]]
      if (.Internal(match(i, currentNetwork[[1]][[j]], F, NULL))){
        # REMOVE i --> j

        # this removes ONLY the first instance of i.
        # but there should not be more than 1 instance
        proposalNetwork[[1]][[j]] <- proposalNetwork[[1]][[j]][-match(i,
                                                  proposalNetwork[[1]][[j]])]

        # update the routes matrix for the proposal
        proposalNetwork[[2]] <- routesRemoveEdge(proposalNetwork[[2]], i, j)

        pr <- prior(proposalNetwork[[1]])
        proposalNetwork[[3]] <- log(pr)
        proposalNetwork[[4]][i, j] <- 0
        proposalNetwork[[5]] <- logNumMHNeighbours(proposalNetwork[[2]],
                                                   proposalNetwork[[4]],
                                                   constraintT)
        if (pr > 0){
          logAccProb <- logScoreOnlineFUN(
                          currentBN          = currentNetwork[[1]],
                          proposalBN         = proposalNetwork[[1]],
                          heads              = j,
                          logScoreParameters = logScoreParameters,
                          cache              = cache,
                          checkInput         = F) +
                        proposalNetwork[[3]] -
                        proposalNetwork[[5]] -
                        currentNetwork[[3]] +
                        currentNetwork[[5]]
        }
        else {
          logAccProb <- -Inf
        }
      }
      else {
        # ADD i --> j

        # this fast unique.default(c(net[[j]], i))
        # the sorting is REQUIRED to canonicalise the network IDs
        # the sorting is fast sort.int(x, method = "quick")
        proposalNetwork[[1]][[j]] <- .Internal(sort(c(i,
                                proposalNetwork[[1]][[j]]), F))

        # update the routes matrix for the proposal
        proposalNetwork[[2]] <- routesAddEdge(proposalNetwork[[2]], i, j)
        pr <- prior(proposalNetwork[[1]])
        proposalNetwork[[3]] <- log(pr)
        proposalNetwork[[4]][i, j] <- 1
        proposalNetwork[[5]] <- logNumMHNeighbours(proposalNetwork[[2]],
                                                   proposalNetwork[[4]],
                                                   constraintT)

        if (pr > 0){
          logAccProb <- logScoreOnlineFUN(
                          currentBN          = currentNetwork[[1]],
                          proposalBN         = proposalNetwork[[1]],
                          heads              = j,
                          logScoreParameters = logScoreParameters,
                          cache              = cache,
                          checkInput         = F) +
                        proposalNetwork[[3]] -
                        proposalNetwork[[5]] -
                        currentNetwork[[3]] +
                        currentNetwork[[5]]
        }
        else {
          logAccProb <- -Inf
        }
      }
    }
    else {
      ## FLIP i --> j

      i <- nonCycleInducingFlips[[select - nNonCycleInducing, 1]]
      j <- nonCycleInducingFlips[[select - nNonCycleInducing, 2]]

      # remove i --> j
      proposalNetwork[[1]][[j]] <- proposalNetwork[[1]][[j]][-match(i,
                                                proposalNetwork[[1]][[j]])]
      proposalNetwork[[4]][i, j] <- 0

      # add j --> i
      proposalNetwork[[1]][[i]] <- .Internal(
                                     sort(c(j, proposalNetwork[[1]][[i]]),
                                          F))
      proposalNetwork[[4]][j, i] <- 1

      proposalNetwork[[2]] <- routesRemoveEdge(proposalNetwork[[2]], i, j)
      proposalNetwork[[2]] <- routesAddEdge(proposalNetwork[[2]], j, i)

      pr <- prior(proposalNetwork[[1]])
      proposalNetwork[[3]] <- log(pr)
      proposalNetwork[[5]] <- logNumMHNeighbours(proposalNetwork[[2]],
                                                 proposalNetwork[[4]],
                                                 constraintT)

      if (pr > 0){
        logAccProb <- logScoreOnlineFUN(
                        currentBN          = currentNetwork[[1]],
                        proposalBN         = proposalNetwork[[1]],
                        heads              = c(j, i),
                        cache              = cache,
                        logScoreParameters = logScoreParameters,
                        checkInput         = F) +
                      proposalNetwork[[3]] -
                      proposalNetwork[[5]] -
                      currentNetwork[[3]] +
                      currentNetwork[[5]]
        }
        else {
          logAccProb <- -Inf
        }
    }

    if (isTRUE(debugAcceptance)) browser()

    if (logAccProb >= 0 || logp < logAccProb){

      # if ACCEPTING the proposal
      currentNetwork <<- proposalNetwork
      nAccepted <<- nAccepted + 1
      if (isTRUE(keepTape)) updateTape(nSteps,
                                       movetype = "mh",
                                       logAccProb,
                                       currentNetwork,
                                       accepted = T)
    }
    else {
      # if REJECTING the proposal
      if (isTRUE(keepTape)) updateTape(nSteps,
                                       movetype = "mh",
                                       logAccProb,
                                       currentNetwork,
                                       accepted = F)
    }
    # return
    # either the logScore of the network
    # or the network
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
