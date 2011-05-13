# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Create a MCMC sampler for Bayesian Networks. The sampler samples Bayesian
#' Networks (ie models).
#'
#' .....
#'
#' @param data The data.
#' @param initial An object of class 'bn'. The starting value of the MCMC.
#' @param prior A function that returns the prior score of the supplied
#'   bn.
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
#' @param logScoreParameters ...
#' @param verbose A logical of length 1, indicating whether verbose output
#'   should be printed.
#' @param keepTape A logical of length 1, indicating whether a full log
#'   ('tape') of the MCMC sampler should be kept. Enabling this
#'   option can be very memory-intensive.
#' @return A function, which when called draws the next sample of the MCMC.
#' @export
BNSamplerBigFlips <- function(data,
                      initial,
                      prior,
                      return      = "network",
                      logScoreFUN = defaultLogScoreFUN(),
                      logScoreParameters = list(hyperparameters = "qi"),
                      verbose     = F,
                      keepTape    = F){
  stopifnot(class(data)           ==   "data.frame",
            all(unlist(lapply(data, class)) == "factor"),
            "bn"                  %in% class(initial),
            is.valid(initial),
            ncol(as.matrix(data)) ==   length(initial),
            is.function(prior),
            return                %in% c("network", "contingency"),
            is.logical(keepTape),
            length(keepTape)      ==   1)

  numberOfNodes <- length(initial)
  nodesSeq <- seq_len(numberOfNodes)
  cache <- new.env(hash = T, size = 10000L)
  logScoreOfflineFUN <- logScoreFUN$offline
  logScoreOnlineFUN <- logScoreFUN$online
  prepareDataFUN <- logScoreFUN$prepare

  # Set up for fast computation of logScore
  logScoreParameters <- prepareDataFUN(data, logScoreParameters, checkInput = F)

  # The current MCMC state is stored a list of the form:
  # currentNetwork[[1]] is the bn
  # currentNetwork[[2]] is the routes matrix
  # currentNetwork[[3]] is the log prior
  currentNetwork <- list(5, mode = "list")
  currentNetwork[[1]] <- initial
  currentScore <- logScoreOfflineFUN(x                  = currentNetwork[[1]],
                                     logScoreParameters = logScoreParameters,
                                     cache = cache)
  currentNetwork[[2]] <- routes(currentNetwork[[1]])
  currentNetwork[[3]] <- log(prior(currentNetwork[[1]]))
  if (!is.valid.prior(currentNetwork[[3]])){
    stop("Initial network has prior with 0 probability.")
  }

  # Set up internal counters and logs etc
  nAccepted <- 0
  nSteps <- 0
  nMHProposals <- 0
  nMJProposals <- 0
  nMHAccepted <- 0
  nMJAccepted <- 0
  nCurrentGraphIsAMode <- 0
  if (return == "contingency"){
    count <- new.env()
  }

  if (isTRUE(keepTape)){
    tapeSizeIncrement <- 500000
    tapeColumns <- c("movetype", "logAccProb", "accepted")
    numberTapeColumns <- length(tapeColumns)
    tape <- matrix(nrow = 0, ncol = numberTapeColumns)
    tapeProposals <- character(length = 0)
    colnames(tape) <- tapeColumns
  }

  moves <- which(row(currentNetwork[[2]]) != col(currentNetwork[[2]]),
                 arr.ind = T)

  function(x, verbose = F, returnDiagnostics = F,
           debugAcceptance = F, returnTape = F){

    if (returnDiagnostics == T){
      return(list(nAccepted            = nAccepted,
                  nSteps               = nSteps,
                  acceptanceRate       = nAccepted/nSteps,
                  nMHProposals         = nMHProposals,
                  nMJProposals         = nMJProposals,
                  nMHAccepted          = nMHAccepted,
                  nMJAccepted          = nMJAccepted,
                  MHAcceptanceRate     = nMHAccepted/nMHProposals,
                  MJAcceptanceRate     = nMJAccepted/nMJProposals,
                  nCurrentGraphIsAMode = nCurrentGraphIsAMode))
    }
    if (returnTape == T){
      if (isTRUE(keepTape)){
        return(data.frame(tape[seq_len(nSteps), ],
                          proposals = tapeProposals[seq_len(nSteps)]))
      }
      else {
        stop("Tape not kept for this MCMC run")
      }
    }
    if (debugAcceptance == T) browser()

    if (isTRUE(keepTape)){
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

    nSteps <<- nSteps + 1

    proposalNetwork <- currentNetwork

    # function for generating a proposal
    # returns the acceptance probability

    movetype <- "mh"
    logp <- log(runif(1, min = 0, max = 1))

    select <- sample.int(nrow(moves), size = 1, replace = T)
    i <- moves[[select, 1]]
    j <- moves[[select, 2]]
    type <- NULL
    if (currentNetwork[[2]][j, i] == 0){
      # the condition is a marginally faster i %in% currentNetwork[[1]][[j]]
      if (.Internal(match(i, currentNetwork[[1]][[j]], F, NULL))){
        # REMOVE i --> j
        type <- 1

        # this removes ONLY the first instance of i.
        # but there should not be more than 1 instance
        proposalNetwork[[1]][[j]] <- proposalNetwork[[1]][[j]][
                                                  -.Internal(match(i,
                                                    proposalNetwork[[1]][[j]],
                                                    F, NULL))]

        # get the prior for the proposal
        pr <- prior(proposalNetwork[[1]])
        proposalNetwork[[3]] <- log(pr)

        if (pr > 0){
          logAccProb <- logScoreOnlineFUN(
                          currentBN          = currentNetwork[[1]],
                          proposalBN         = proposalNetwork[[1]],
                          heads              = j,
                          logScoreParameters = logScoreParameters,
                          cache              = cache,
                          checkInput         = F) +
                        proposalNetwork[[3]] -
                        currentNetwork[[3]]
        }
        else {
          logAccProb <- -Inf
        }
      }
      else {
        # ADD i --> j
        type <- 2

        # this fast unique.default(c(net[[j]], i))
        # the sorting is REQUIRED to canonicalise the network IDs
        # the sorting is fast sort.int(x, method = "quick")

        proposalNetwork[[1]][[j]] <- .Internal(sort(c(i,
                                           proposalNetwork[[1]][[j]]), F))

        # get the prior for the proposal
        pr <- prior(proposalNetwork[[1]])
        proposalNetwork[[3]] <- log(pr)

        if (pr > 0){
          logAccProb <- logScoreOnlineFUN(
                          currentBN          = currentNetwork[[1]],
                          proposalBN         = proposalNetwork[[1]],
                          heads              = j,
                          logScoreParameters = logScoreParameters,
                          cache              = cache,
                          checkInput         = F) +
                        proposalNetwork[[3]] -
                        currentNetwork[[3]]
        }
        else {
          logAccProb <- -Inf
        }
      }
    }
    else {
      ## FLIP i --> j
      # or sometimes a bigger flip
      type <- 3

      onRoute <- nodesSeq[currentNetwork[[2]][,i] & currentNetwork[[2]][j,]]

      # loop over each node on the path
      edgesOnPath <- unlist(lapply(onRoute, function(node) {
        parentsOnRoute <- .Internal(unique(onRoute[.Internal(
          match(proposalNetwork[[1]][[node]], onRoute, F, NULL))],
          F, F))
        lapply(parentsOnRoute, function(parent){
          c(node, parent)
        })
      }), recursive = F)

      for (rowNum in seq_len(length(edgesOnPath))){
        pair <- edgesOnPath[[rowNum]]
        proposalNetwork[[1]][[pair[1]]] <- proposalNetwork[[1]][[pair[1]]][
                                        -.Internal(match(pair[2],
                                        proposalNetwork[[1]][[pair[1]]],
                                        F, F))]
        proposalNetwork[[1]][[pair[2]]] <- .Internal(sort(c(pair[1],
                                     proposalNetwork[[1]][[pair[2]]]), F))
      }

      # get the prior for the proposal
      pr <- prior(proposalNetwork[[1]])
      proposalNetwork[[3]] <- log(pr)

      if (pr > 0){
        logAccProb <- logScoreOnlineFUN(
                        currentBN          = currentNetwork[[1]],
                        proposalBN         = proposalNetwork[[1]],
                        heads              = onRoute,
                        cache              = cache,
                        logScoreParameters = logScoreParameters,
                        checkInput         = F) +
                      proposalNetwork[[3]] -
                      currentNetwork[[3]]
        }
        else {
          logAccProb <- -Inf
        }
    }

    if (isTRUE(keepTape)){
      tape[nSteps, 2] <<- logAccProb
      tapeProposals[nSteps] <<- as.character(proposalNetwork[[1]], pretty = T)
    }

    if (debugAcceptance == T) browser()

    if (logAccProb >= 0 || logp < logAccProb){
      # if ACCEPTING the proposal

      # update the routes matrix
      if (type == 1){
        proposalNetwork[[2]] <- proposalNetwork[[2]] -
                                  proposalNetwork[[2]][, i] %*%
                                  .Internal(t.default(
                                    proposalNetwork[[2]][j, ]))
      }
      else if (type == 2){
        proposalNetwork[[2]] <- proposalNetwork[[2]] +
                                  proposalNetwork[[2]][, i] %*%
                                  .Internal(t.default(
                                    proposalNetwork[[2]][j, ]))
      }
      else if (type == 3){
        # can't do this in one because the intermediate graphs
        # may well be cyclic!
        for (rowNum in seq_len(length(edgesOnPath))){
          pair <- edgesOnPath[[rowNum]]
          proposalNetwork[[2]] <- proposalNetwork[[2]] -
                                    proposalNetwork[[2]][, pair[2]] %*%
                                      .Internal(t.default(
                                        proposalNetwork[[2]][pair[1], ]))
        }

        for (rowNum in seq_len(length(edgesOnPath))){
          pair <- edgesOnPath[[rowNum]]
          proposalNetwork[[2]] <- proposalNetwork[[2]] +
                                    proposalNetwork[[2]][, pair[1]] %*%
                                    .Internal(t.default(
                                      proposalNetwork[[2]][pair[2], ]))
        }
      }

      currentNetwork <<- proposalNetwork
      nAccepted <<- nAccepted + 1

      # return
      # either the logScore of the network
      # or the network
      currentNetwork[[1]]
    }
    else {
      # if REJECTING the proposal

      # return
      # either the logScore of the network
      # or the network
      currentNetwork[[1]]
    }
  }
}