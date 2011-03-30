# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Identify neighbouring graphs
#'
#' Finds which of the supplied graphs are not neighbours. Neighbouring
#' graphs are removed according to a heuristic involving scores.
#'
#' @param x A list of graphs
#' @param logScores A numeric vector containing the log score of the 
#'   corresponding graph in \code{x}
#' @param verbose Should verbose output be shown. A logical of length 1.
#' @param head The rank of the worst graph to return.
#' @return A numeric vector, giving the indicies \code{x} which are 
#'   not neighbours.
#' @export
whichGraphsNotNeighbours <- function(x, logScores, verbose = F, head){
  stopifnot("parental.list"     %in% class(x),
            class(logScores)  == "numeric",
            length(logScores) == length(x),
            class(verbose)      == "logical",
            length(verbose)     == 1)
  included <- seq_along(x)
  numberOfModes <- length(x)
  modesDifference <- matrix(NA, nrow = numberOfModes, ncol = numberOfModes)
  modesDifference[upper.tri(modesDifference)] <- -1

  # unique() works on lists if everything is identical
  # e.g.
  #     duplicated(list(list(c(1,2),3L), list(4,5), list(c(1,2),3L)))
  # but..
  #     duplicated(list(list(c(1,2),3L), list(4,5), list(c(1,2),3)))
  # and..
  #     duplicated(list(list(c(2,1),3L), list(4,5), list(c(1,2),3)))
  #

  # but here we should be OK if all the parentals in x are valid!
  whichDuplicated <- which(duplicated(x))
  included <- setdiff(included, whichDuplicated)

  select <- matrix(F, nrow = numberOfModes, ncol = numberOfModes)
  select[whichDuplicated, ] <- T
  select[, whichDuplicated] <- T
  modesDifference[select & upper.tri(select)] <- 0

  if (!missing(head)){
    whichNotInHead <- which(rank(-logScores, ties.method = "min") > head)
    included <- setdiff(included, whichNotInHead)

    select <- matrix(F, nrow = numberOfModes, ncol = numberOfModes)
    select[whichNotInHead, ] <- T
    select[, whichNotInHead] <- T
    modesDifference[select & upper.tri(select)] <- 0
  }

  whichToCheck <- which(modesDifference == -1, arr.ind = T)

  if (verbose){
    cat("Computing the number of moves between the modes...\n")
    cat(nrow(whichToCheck), "operations")
  }
  for (row in seq_len(nrow(whichToCheck))){
    i <- whichToCheck[row, 1]
    j <- whichToCheck[row, 2]

    modesDifference[i, j] <- numberOfMovesBetweenIgnoringCycles(
                              x[[i]], x[[j]], allowFlips = T)
  }

  if (verbose){
    cat("Removing adjacent modes...\n")
  }
  while (sum(modesDifference == 1, na.rm = T) != 0){
    whichModesAreNeighbours <- which(modesDifference == 1, arr.ind = T)

    whichModesAreNeighboursNumeric <- as.numeric(whichModesAreNeighbours)

    # there is some scope for some kind of heuristic that says that

    # if removing just one mode would make everything OK, and that mode is
    # not terribly great anyway then we just remove it.
    #if (any(table(whichModesAreNeighboursNumeric) == length(included) - 1)){
    #}

    modesThatAreAdjacentToSomething <- unique(whichModesAreNeighboursNumeric)

    modeToRemove <- modesThatAreAdjacentToSomething[which(order(
      logScores[modesThatAreAdjacentToSomething]) == 1)]

    included <- setdiff(included, modeToRemove)
    modesDifference[modeToRemove, ] <- 0
    modesDifference[, modeToRemove] <- 0
  }

  if (verbose){
    cat(paste(length(included), "remain after removing neighbours\n"))
  }
  included
}

#' Mode-jumping MCMC sampler for Bayesian Networks.
#' 
#' The sampler samples Bayesian Networks (ie models).
#'
#' @param data The data.
#' @param initial An object of class 'bn'. The starting value of the MCMC.
#' @param prior A function that returns the prior score of the supplied bn.
#' @param return Either "network" or "contingency".
#' @param logScoreFUN A list of three elements:
#'   \describe{
#'     \item{offline}{A function that computes the logScore of a Bayesian 
#'                    Network}
#'     \item{online}{A function that incrementally computes the logScore of a 
#'                   Bayesian Network}
#'     \item{prepare}{A function that prepares the data, and any further 
#'                    pre-computation required by the logScore functions.}
#'   }
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
#' @param modejumping Either a logical of length 1, or a list. When no
#'   mode-jumping is desired, use modejumping = F. For mode-jumping, use a 
#'   list with the following components:
#'   \describe{
#'     \item{modes}{A bn.list of modes.}
#'     \item{modesLogScores}{Optionally, a numeric vector, containing the 
#'                           logScores of the modes.}
#'     \item{checkModesAcyclic}{A logical of length 1}
#'     \item{modesPreFiltered}{A logical of length 1}
#'     \item{modeJumpingProbability}{A numeric of length 1. Default 0.25}
#'     \item{dontCheckModesValid}{A logical of length 1}
#'   }
#' @param verbose A logical of length 1, indicating whether verbose
#'  output should be printed.
#' @param keepTape A logical of length 1, indicating whether a full log
#'   (\code{tape}) of the MCMC sampler should be kept. Enabling this option 
#'   can be very memory-intensive.
#' @return A function, which when called draws the next sample of the MCMC.
#' @export
BNSamplerMJ <- function(data,
                        initial,
                        prior,
                        return             = "network",
                        logScoreFUN        = defaultLogScoreFUN(),
                        logScoreParameters = list(hyperparameters = "qi"),
                        constraint         = NULL,
                        modejumping        = F,
                        verbose            = F,
                        keepTape           = F){
  stopifnot("bn"                  %in% class(initial),
            is.valid(initial),
            ncol(as.matrix(data)) ==   length(initial),
            is.function(prior),
            return                %in% c("network", "contingency", "sha256"),
            is.list(modejumping) || is.logical(modejumping),
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
  usingConstraint <- any(notdiag(constraint) > 0)
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
  nMJProposals <- 0
  nMHAccepted <- 0
  nMJAccepted <- 0
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

  # Set up mode-jumping
  if (is.list(modejumping)){
    modes <- modejumping$modes
    modesLogScores <- modejumping$modesLogScores
    modeJumpingProbability <- modejumping$modeJumpingProbability
    dontCheckModesValid <- modejumping$dontCheckModesValid
    dontCheckModesConstraint <- modejumping$dontCheckModesConstraint
    if (is.null(modeJumpingProbability)) modeJumpingProbability <- 0.25
    checkModesAcyclic <- modejumping$checkModesAcyclic
    if (is.null(checkModesAcyclic)) checkModesAcyclic <- T
    modesPreFiltered <- modejumping$modesPreFiltered
    if (is.null(modesPreFiltered)) modesPreFiltered <- F
    if (is.null(dontCheckModesValid)) dontCheckModesValid <- F
    if (is.null(dontCheckModesConstraint)) dontCheckModesConstraint <- F

    if (!dontCheckModesValid){
      stopifnot(all(sapply(modes, is.valid) == T))
    }

    if (is.null(modesLogScores)){
      modesLogScores <- sapply(modes, function(mode){
        logScoreOfflineFUN(x                  = mode,
                           logScoreParameters = logScoreParameters,
                           cache              = cache)
      })
    }

    if (!dontCheckModesConstraint){
      if (usingConstraint){
        if (!any(unlist(lapply(modes, satisfiesConstraint, constraint = constraint)))){
          stop("At least one of the modes does not satisfy the constraint")
        }
      }
    }

    if (isTRUE(checkModesAcyclic)){
      if (verbose){
        cat("Checking modes for cycles\n")
      }
      # check all the modes are acyclic
      if (!any(unlist(lapply(modes, checkAcyclic)))){
        stop("At least one of the modes is not acyclic")
      }
    }

    numberOfModes <- length(modes)
    modesSeq <- seq_along(modes)

    # remove neighbouring modes
    if (!isTRUE(modesPreFiltered)){
      if (verbose){
        cat("Computing number of moves between modes\n")
      }
      included <- whichGraphsNotNeighbours(modes, modesLogScores, verbose)
      modes <- modes[included]
      modesLogScores <- modesLogScores[included]
    }
    numberOfModes <- length(modes)

    if (numberOfModes == 1){
      stop("Only one non-adjacent, unique mode.")
    } else if (numberOfModes <= 10){
      warning("Number of non-adjacenct, unique modes (there are ",
                          numberOfModes, ") included is low.")
    }

    modesID <- sapply(modes, fastid)
    currentGraphIsAMode <- fastid(currentNetwork[[1]]) %in% modesID
  }

  returnDiagnostics <- function(){
    list(nAccepted            = nAccepted,
         nSteps               = nSteps,
         acceptanceRate       = nAccepted/nSteps,
         nMHProposals         = nMHProposals,
         nMJProposals         = nMJProposals,
         nMHAccepted          = nMHAccepted,
         nMJAccepted          = nMJAccepted,
         MHAcceptanceRate     = nMHAccepted/nMHProposals,
         MJAcceptanceRate     = nMJAccepted/nMJProposals,
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
                         currentGraphIsAMode,
                         proposalGraphIsAMode,
                         accepted){
    tape[nSteps, 1] <<- if (movetype == "mh") 0 else 1
    tape[nSteps, 2] <<- logAccProb
    tapeProposals[nSteps] <<- as.character(currentNetwork[[1]], pretty = T)
    if (is.list(modejumping)){
      tape[nSteps, 4] <<- currentGraphIsAMode
      tape[nSteps, 5] <<- proposalGraphIsAMode
    }
  }

  function(x, verbose = F, returnDiagnostics = F,
           debugAcceptance = F, returnTape = F, burnin = 0){
    if (isTRUE(returnDiagnostics)) return(returnDiagnostics())
    if (isTRUE(returnTape)) return(returnTape())
    if (isTRUE(debugAcceptance)) browser()
    if (isTRUE(keepTape)) lengthenTape()

    nSteps <<- nSteps + 1
    proposalNetwork <- currentNetwork

    # function for generating a proposal
    # returns the acceptance probability

    movetype <- "mh"
    logp <- log(runif(1, min = 0, max = 1))

    if (is.list(modejumping)){
      proposalGraphIsAMode <- F
      if (isTRUE(currentGraphIsAMode)){
        nCurrentGraphIsAMode <<- nCurrentGraphIsAMode + 1

        if (runif(1, min = 0, max = 1) < modeJumpingProbability){
          # propose mode-jumping move
          movetype <- "mj"
          nMJProposals <<- nMJProposals + 1
          proposalGraphIsAMode <- T

          currentID <- fastid(currentNetwork[[1]])
          whichMode <- match(currentID, modesID)

          #otherModes <- modes[-whichMode]
          #proposalNetwork[[1]] <- otherModes[[sample(seq_along(otherModes), 1)]]

          nModes <- length(modes)
          modesSeqMinus1 <- seq_len(nModes - 1)
          samp <- sample(modesSeqMinus1, 1)
          if (samp >= whichMode){
            samp <- samp + 1
          }
          proposalNetwork[[1]] <- modes[[samp]]

          proposalNetwork[[2]] <- routes(proposalNetwork[[1]])
          pr <- prior(proposalNetwork[[1]])
          proposalNetwork[[3]] <- log(pr)
          proposalNetwork[[4]] <- as.adjacency(proposalNetwork[[1]])
          proposalNetwork[[5]] <- logNumMHNeighbours(proposalNetwork[[2]],
                                                     proposalNetwork[[4]],
                                                     constraintT)

          logScoreOfflineFUN(x                  = proposalNetwork[[1]],
                             logScoreParameters = logScoreParameters,
                             cache              = cache,
                             checkInput         = F)

          if (pr > 0){
            # don't use the neighbourhood size here
            logAccProb <- logScoreOnlineFUN(
                            currentBN          = currentNetwork[[1]],
                            proposalBN         = proposalNetwork[[1]],
                            heads              = nodesSeq,
                            logScoreParameters = logScoreParameters,
                            cache              = cache,
                            checkInput         = F) +
                          proposalNetwork[[3]] - currentNetwork[[3]]
          }
          else {
            logAccProb <- -Inf
          }
        }
      }
    }
    if (movetype == "mh"){
      # propose mh move
      nMHProposals <<- nMHProposals + 1

      # count the number of proposals and select one
      canAddOrRemove <- currentNetwork[[2]] == 0 & constraintT == 0
      canFlip <- currentNetwork[[2]] == 1 &
                 currentNetwork[[4]] == 1 &

                 constraintT == 0

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

      if (is.list(modejumping)){
        proposalID <- fastid(proposalNetwork[[1]])
        if (proposalID %in% modesID){
          proposalGraphIsAMode <- T
          if (runif(1, min = 0, max = 1) < modeJumpingProbability){
            # REJECT the M-H proposal
            logp <- 1
            logAccProb <- -1
          }
        }
      }
    }

    if (isTRUE(debugAcceptance)) browser()

    if (logAccProb >= 0 || logp < logAccProb){

      # if ACCEPTING the proposal
      currentNetwork <<- proposalNetwork
      nAccepted <<- nAccepted + 1

      if (isTRUE(keepTape)) updateTape(nSteps,
                                       currentNetwork,
                                       movetype,
                                       logAccProb,
                                       currentGraphIsAMode,
                                       proposalGraphIsAMode,
                                       accepted = T)

      if (is.list(modejumping)){
        if (movetype == "mh"){
          nMHAccepted <<- nMHAccepted + 1
        }
        else {
          nMJAccepted <<- nMJAccepted + 1
        }

        if (proposalGraphIsAMode){
          currentGraphIsAMode <<- T
        }
        else {
          currentGraphIsAMode <<- F
        }
      }
    }
    else {
      # if REJECTING the proposal
      if (isTRUE(keepTape)) updateTape(nSteps,
                                       currentNetwork,
                                       movetype,
                                       logAccProb,
                                       currentGraphIsAMode,
                                       proposalGraphIsAMode,
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