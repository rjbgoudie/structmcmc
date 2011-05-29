# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Undocumented.
#'
#' method description
#'
#' @param samplers ...
#' @param data ...
#' @param n ...
#' @param temperatures ...
#' @param pswap ...
#' @param burnin ...
#' @param verbose ...
#' @param currentNetwork BUG found by R-check
#' @export
#' @seealso \code{\link{BNSamplerPT}}, \code{\link{draw}}
drawPT <- function(samplers,
                   data,
                   n,
                   temperatures,
                   pswap = 0.2,
                   burnin = 0,
                   verbose = T,
                   currentNetwork){

  nTemperatures <- length(temperatures)
  samplersSeq <- seq_along(samplers)
  # Note: the initial graph is NOT returned at the moment
  samplesl <- lapply(samplersSeq, function(i){
    vector("list", n)
  })
  if (verbose){
    cat("Drawing", n, "samples. Counting 10000s: ")
  }
  nStandardMoves <- 0
  nTemperatureFlipMoves <- 0
  nTemperatureFlipMovesAccepted <- 0

  for (i in seq_len(n)){
    if (verbose && i %% 100 == 0){
      cat(i/100 , ", ", sep = "")
    }

    u <- runif(n = 1, min = 0, max = 1)

    if (u < 1 - pswap){
      nStandardMoves <- nStandardMoves + 1
      for (s in samplersSeq){
        out <- samplers[[s]](i)
        class(out) <- c("bn", "parental")
        samplesl[[s]][[i]] <- out
      }
    } else {
      nTemperatureFlipMoves <- nTemperatureFlipMoves + 1

      a <- runif(n = 1, min = 1, max = nTemperatures)
      chaini <- floor(a)
      chainj <- ceiling(a)

      prevj <- evalq(currentNetwork[[1]], envir = environment(samplers[[chainj]]))
      previ <- evalq(currentNetwork[[1]], envir = environment(samplers[[chaini]]))

      rho <- (temperatures[chaini] * logScoreMultDir(prevj, data) +
      temperatures[chainj] * logScoreMultDir(previ, data)) -
      (temperatures[chaini] * logScoreMultDir(previ, data) +
      temperatures[chainj] * logScoreMultDir(prevj, data))

      v <- log(runif(n = 1, min = 0, max = 1))

      if (rho >= 0 || v < rho){
        nTemperatureFlipMovesAccepted <- nTemperatureFlipMovesAccepted + 1
        samplesl[[chaini]][[i]] <- prevj
        samplesl[[chainj]][[i]] <- previ

        for (s in setdiff(samplersSeq, c(chaini, chainj))){
          samplesl[[s]][[i]] <- evalq(currentNetwork[[1]], envir = environment(samplers[[s]]))
        }
      } else {
        for (s in samplersSeq){
          samplesl[[s]][[i]] <- evalq(currentNetwork[[1]], envir = environment(samplers[[s]]))
        }
      }
    }
  }

  if (verbose){
    cat("\n")
    cat("Normal moves:", nStandardMoves,
        "\nTempFlipMoves:", nTemperatureFlipMoves,
        "\nTempFlipMovesAccepted", nTemperatureFlipMovesAccepted)
  }

  samplesl <- lapply(samplesl, function(samples){
    class(samples) <- c("mcmcbn", "bn.list", "parental.list")
    samples
  })
  samplesl[[1]]
}

#' Undocumented.
#'
#' method description
#'
#' @param data ...
#' @param initial ...
#' @param prior ...
#' @param return ...
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
#' @param logScoreParameters ...
#' @param constraint ...
#' @param modejumping ...
#' @param verbose ...
#' @param keepTape ...
#' @param temp ...
#' @param cache ...
#' @export
#' @seealso \code{\link{drawPT}}. \code{\link{BNSampler}},
#'   \code{\link{BNSamplerMJ}}, \code{\link{BNGibbsSampler}},
#'   \code{\link{BNSamplerBigFlips}}, \code{\link{BNSamplerGrzeg}}
BNSamplerPT <- function(data,
                      initial,
                      prior,
                      return      = "network",
                      logScoreFUN = logScoreMultDirFUN(),
                      logScoreParameters = list(hyperparameters = "qi"),
                      constraint  = NULL,
                      modejumping = F,
                      verbose     = F,
                      keepTape    = F,
                      temp        = 1,
                      cache       = new.env(hash = T, size = 10000L)){
  stopifnot(class(data)           ==   "data.frame",
            all(unlist(lapply(data, class)) == "factor"),
            "bn"                  %in% class(initial),
            is.valid(initial),
            ncol(as.matrix(data)) ==   length(initial),
            is.function(prior),
            return                %in% c("network", "contingency"),
            is.list(modejumping) || is.logical(modejumping),
            is.logical(keepTape),
            length(keepTape)      ==   1)

  numberOfNodes <- length(initial)
  nodesSeq <- seq_len(numberOfNodes)

  logScoreOfflineFUN <- logScoreFUN$offline
  logScoreOnlineFUN <- logScoreFUN$online
  prepareDataFUN <- logScoreFUN$prepare

  # constraints
  if (is.null(constraint)){
    usingConstraint <- F
    constraint <- matrix(0, numberOfNodes, numberOfNodes)
    constraintT <- t(constraint)
  }
  else {
    if (is.list(modejumping)){
      warning("Mode jumping may not work with constraints at the moment.")
    }
    stopifnot(class(constraint)    ==   "matrix",
              all(constraint       %in% c(-1, 0, 1)),
              all(diag(constraint) ==   0))
    usingConstraint <- T

    # check initial meet constraint
    if (!satisfiesConstraint(initial, constraint)){
      stop("Initial network does not satisfy constraint")
    }
    constraintT <- t(constraint)
  }

  # Set up for fast computation of logScoreMultDir
  logScoreParameters <- prepareDataFUN(data,
                                       logScoreParameters,
                                       checkInput = F)

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

  # Set up mode-jumping
  if (is.list(modejumping)){
    modes <- modejumping$modes
    modesLogScores <- modejumping$modesLogScores
    modeJumpingProbability <- modejumping$modeJumpingProbability
    if (is.null(modeJumpingProbability)) modeJumpingProbability <- 0.25
    checkModesAcyclic <- modejumping$checkModesAcyclic
    if (is.null(checkModesAcyclic)) checkModesAcyclic <- T
    modesPreFiltered <- modejumping$modesPreFiltered
    if (is.null(modesPreFiltered)) modesPreFiltered <- F

    stopifnot(all(sapply(modes, is.valid) == T))

    if (is.null(modesLogScores)){
      modesLogScores <- sapply(modes, function(mode){
        logScoreOfflineFUN(x                  = mode,
                           logScoreParameters = logScoreParameters,
                           cache              = cache)
      })
    }

    if (usingConstraint){
      if (!any(unlist(lapply(modes, satisfiesConstraint, constraint = constraint)))){
        stop("At least one of the modes does not satisfy the constraint")
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

    modesID <- lapply(modes, fastid)
    currentGraphIsAMode <- fastid(currentNetwork[[1]]) %in% modesID
  }

  sampler <- function(x, verbose = F, returnDiagnostics = F,
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
          whichMode <- which(modesID == currentID)
          otherModes <- modes[-whichMode]
          proposalNetwork[[1]] <- otherModes[[sample(seq_along(otherModes), 1)]]
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
            logAccProb <- temp * (logScoreOnlineFUN(
                            currentBN          = currentNetwork[[1]],
                            proposalBN         = proposalNetwork[[1]],
                            heads              = nodesSeq,
                            logScoreParameters = logScoreParameters,
                            cache              = cache,
                            checkInput         = F) +
                          proposalNetwork[[3]] - currentNetwork[[3]])
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
            logAccProb <- temp * (logScoreOnlineFUN(
                            currentBN          = currentNetwork[[1]],
                            proposalBN         = proposalNetwork[[1]],
                            heads              = j,
                            logScoreParameters = logScoreParameters,
                            cache              = cache,
                            checkInput         = F) +
                          proposalNetwork[[3]] -
                          proposalNetwork[[5]] -
                          currentNetwork[[3]] +
                          currentNetwork[[5]])
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
            logAccProb <- temp * (logScoreOnlineFUN(
                            currentBN          = currentNetwork[[1]],
                            proposalBN         = proposalNetwork[[1]],
                            heads              = j,
                            logScoreParameters = logScoreParameters,
                            cache              = cache,
                            checkInput         = F) +
                          proposalNetwork[[3]] -
                          proposalNetwork[[5]] -
                          currentNetwork[[3]] +
                          currentNetwork[[5]])
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
          logAccProb <- temp * (logScoreOnlineFUN(
                          currentBN          = currentNetwork[[1]],
                          proposalBN         = proposalNetwork[[1]],
                          heads              = c(j, i),
                          cache              = cache,
                          logScoreParameters = logScoreParameters,
                          checkInput         = F) +
                        proposalNetwork[[3]] -
                        proposalNetwork[[5]] -
                        currentNetwork[[3]] +
                        currentNetwork[[5]])
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

    if (isTRUE(keepTape)){
      tape[nSteps, 1] <<- if (movetype == "mh") 0 else 1
      tape[nSteps, 2] <<- logAccProb
      tapeProposals[nSteps] <<- as.character(proposalNetwork[[1]], pretty = T)
    }

    if (debugAcceptance == T) browser()

    if (logAccProb >= 0 || logp < logAccProb){

      # if ACCEPTING the proposal
      currentNetwork <<- proposalNetwork
      nAccepted <<- nAccepted + 1

      if (isTRUE(keepTape)){
        tape[nSteps, 3] <<- 1
      }

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

      # return
      # either the logScore of the network
      # or the network
      if (return == "network"){
        currentNetwork[[1]]
      }
      else {
        id <- do.call("paste", list(currentNetwork[[1]],
                              sep = "", collapse = ","))
        if (is.null(count[[id]])){
          count[[id]] <- 0
        }
        else {
          count[[id]] <- count[[id]] + 1
        }
        NULL
      }
    }
    else {
      # if REJECTING the proposal

      if (isTRUE(keepTape)){
        tape[nSteps, 3] <<- 0
      }

      # return
      # either the logScore of the network
      # or the network
      if (return == "network"){
        currentNetwork[[1]]
      }
      else {
        id <- do.call("paste", list(currentNetwork[[1]],
                                  sep = "", collapse = ","))
        if (is.null(count[[id]])){
          count[[id]] <- 0
        }
        else {
          count[[id]] <- count[[id]] + 1
        }
        NULL
      }
    } 
  }
  class(sampler) <- c("sampler", "function")
  sampler
}
