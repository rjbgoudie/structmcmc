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
#' @param possibleParents ...
#' @param head ...
#' @param proposal ...
#' @param current ...
#' @param logScoreOnlineFUN ...
#' @param logScoreParameters ...
#' @param cache ...
#' @export
scoreParents <- function(possibleParents,
                         head,
                         proposal,
                         current,
                         logScoreOnlineFUN,
                         logScoreParameters,
                         cache){
  sapply(possibleParents, function(parents){
    new <- proposal
    new[[head]] <- parents
    logScoreOnlineFUN(
                    currentBN          = current,
                    proposalBN         = new,
                    heads              = head,
                    logScoreParameters = logScoreParameters,
                    cache              = cache,
                    checkInput         = F)
  })
}

#' Undocumented.
#'
#' method description
#'
#' @param currentNetwork ...
#' @param proposalNetwork ...
#' @param node ...
#' @param required ...
#' @param logScoreOnlineFUN ...
#' @param logScoreParameters ...
#' @param cache ...
#' @export
localPartitionFunction <- function(currentNetwork,
                                   proposalNetwork,
                                   node,
                                   required = integer(0),
                                   logScoreOnlineFUN,
                                   logScoreParameters,
                                   cache){
  nonDescendants <- which(proposalNetwork[[2]][node, ] == 0)
  # require the previous tail to be in the parent set
  nonDescendants <- setdiff(nonDescendants, required)

  possibleParents <- enumerateParents(nonDescendants,
                                       required = required)
  possibleParentsScores <- scoreParents(possibleParents = possibleParents,
                                         head = node,
                                         proposal = proposalNetwork[[1]],
                                         current = currentNetwork[[1]],
                                         logScoreOnlineFUN = logScoreOnlineFUN,
                                         logScoreParameters = logScoreParameters,
                                         cache = cache)

  partition <- logsumexp(possibleParentsScores)
  list(partition = partition,
       weights = exp(possibleParentsScores - partition),
       possibleParents = possibleParents)
}

#' Undocumented.
#'
#' method description. Not tested.
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
#' @param logScoreParameters ...
#' @param constraint ...
#' @param modejumping NOT IMPLEMENTED
#' @param verbose ...
#' @param keepTape ...
#' @param grzeg ...
#' @export
#' @seealso \code{\link{BNSampler}}, \code{\link{BNSamplerPT}},
#'   \code{\link{BNSamplerMJ}}, \code{\link{BNGibbsSampler}}.
#'   Internally uses \code{\link{localPartitionFunction}}, and 
#'   \code{\link{scoreParents}}.
BNSamplerGrzeg <- function(data,
                      initial,
                      prior,
                      return      = "network",
                      logScoreFUN = defaultLogScoreFUN(),
                      logScoreParameters = list(hyperparameters = "qi"),
                      constraint  = NULL,
                      modejumping = F,
                      verbose     = F,
                      keepTape    = F,
                      grzeg       = list(grzegMoveProbability = 0.5)){
  stopifnot("bn"                  %in% class(initial),
            is.valid(initial),
            ncol(as.matrix(data)) ==   length(initial),
            is.function(prior),
            return                %in% c("network", "contingency"),
            is.list(modejumping) || is.logical(modejumping),
            is.logical(keepTape),
            length(keepTape)      ==   1)

  numberOfNodes <- length(initial)
  nodesSeq <- seq_len(numberOfNodes)
  cache <- new.env(hash = T, size = 10000L)
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

  # Set up for fast computation of logScore
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
  currentGraphIsAMode <- F
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

  if (!is.null(grzeg)){
    grzegMoveProbability <- grzeg$grzegMoveProbability
  }

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

    nParents <- sapply(currentNetwork[[1]], length)
    nCurrentEdges <- sum(nParents)

    if (nCurrentEdges > 0){
      # double check this
      if (runif(1, min = 0, max = 1) < grzegMoveProbability){
        movetype <- "grzeg"

        # choose which edge to flip etc
        whichEdge <- sample.int(nCurrentEdges, size = 1)
        head <- match(TRUE, !cumsum(nParents) < whichEdge)
        tail <- unlist(currentNetwork[[1]])[whichEdge]

        # remove the parents of tail
        # ie update proposalNetwork to be
        # M^{\tilde}_{\crosscircle} == M\{parents of tail} by (26)
        for (i in currentNetwork[[1]][[tail]]){
          proposalNetwork[[2]] <- routesRemoveEdge(proposalNetwork[[2]], i, tail)
        }
        proposalNetwork[[1]][[tail]] <- integer(0)
        proposalNetwork[[4]][, tail] <- 0

        # compute log(Z(X_{i} | M^{tilde}_{\crosscircle}))
        p4 <- localPartitionFunction(currentNetwork,
                                     proposalNetwork,
                                     node = tail,
                                     required = integer(0),
                                     logScoreOnlineFUN,
                                     logScoreParameters,
                                     cache)
        logZtail <- p4$partition

        # remove the parents of head
        # ie update proposalNetwork to be
        # M_{\crosscircle} := M^{\tilde}_{\dotcircle}
        for (i in currentNetwork[[1]][[head]]){
          proposalNetwork[[2]] <- routesRemoveEdge(proposalNetwork[[2]], i, head)
        }
        proposalNetwork[[1]][[head]] <- integer(0)
        proposalNetwork[[4]][, head] <- 0

        # compute log(Z*(X_{i} | M_{\dotcircle}, X_{j}))
        p1 <- localPartitionFunction(currentNetwork,
                                     proposalNetwork,
                                     node = tail,
                                     required = head,
                                     logScoreOnlineFUN,
                                     logScoreParameters,
                                     cache)
        logZstarTail <- p1$partition
        weights1 <- p1$weights
        possibleParents1 <- p1$possibleParents

        # sample from Q(\pi^{\tilde}_{tail} | M_{\dotcircle}, X_{j}))
        w1 <- sample.int(length(weights1), size = 1, prob = weights1)
        newParents1 <- sort.int(possibleParents1[[w1]])

        # compute log(Z*(X_{j} | M^{\tilde}_{\dotcircle}, X_{i}))
        p2 <- localPartitionFunction(currentNetwork,
                                     proposalNetwork,
                                     node = head,
                                     required = tail,
                                     logScoreOnlineFUN,
                                     logScoreParameters,
                                     cache)
        logZstarHead <- p2$partition

        # update proposalNetwork to be M_{\crosscircle}
        proposalNetwork[[1]][[tail]] <- newParents1
        proposalNetwork[[4]][newParents1, tail] <- 1
        for (i in newParents1){
          proposalNetwork[[2]] <- routesAddEdge(proposalNetwork[[2]], i, tail)
        }

        # compute log(Z(X_{j} | M_{\crosscircle}))
        p3 <- localPartitionFunction(currentNetwork,
                                     proposalNetwork,
                                     node = head,
                                     required = integer(0),
                                     logScoreOnlineFUN,
                                     logScoreParameters,
                                     cache)
        logZhead <- p3$partition
        weights2 <- p3$weights
        possibleParents2 <- p3$possibleParents

        # sample from Q(\pi^{\tilde} | M_{\crosscircle}))
        w2 <- sample.int(length(weights2), size = 1, prob = weights2)
        newParents2 <- sort.int(possibleParents2[[w2]])

        # update proposalNetwork to be M^{\tilde}
        proposalNetwork[[1]][[head]] <- newParents2
        proposalNetwork[[4]][newParents2, head] <- 1
        for (i in newParents2){
          proposalNetwork[[2]] <- routesAddEdge(proposalNetwork[[2]], i, head)
        }

        # compute log(N^{tilde}^{\cross})
        nProposalEdges <- length(unlist(proposalNetwork[[1]]))

        # acceptance probability is
        # log(N^{\cross}) - log(N^{tilde}^{\cross})
        # log(Z*(X_{i} | M_{\dotcircle}, X_{j})) -
        # log(Z*(X_{j} | M^{\tilde}_{\dotcircle}, X_{i})) +
        # log(Z(X_{j} | M_{\crosscircle})) -
        # log(Z(X_{i} | M^{tilde}_{\crosscircle}))
        logAccProb <- log(nCurrentEdges) - log(nProposalEdges) +
                      logZstarTail - logZstarHead +
                      logZhead - logZtail
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

      # if (is.list(modejumping)){
      #   proposalID <- modeID(proposalNetwork[[1]])
      #   if (proposalID %in% modesID){
      #     proposalGraphIsAMode <- T
      #     if (runif(1, min = 0, max = 1) < grzegMoveProbability){
      #       # REJECT the M-H proposal
      #       logp <- 1
      #       logAccProb <- -1
      #     }
      #   }
      # }
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

      # if (is.list(modejumping)){
      #         if (movetype == "mh"){
      #           nMHAccepted <<- nMHAccepted + 1
      #         }
      #         else {
      #           nMJAccepted <<- nMJAccepted + 1
      #         }
      #        
      #         if (proposalGraphIsAMode){
      #           currentGraphIsAMode <<- T
      #         }
      #         else {
      #           currentGraphIsAMode <<- F
      #         }
      #       }

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
}
