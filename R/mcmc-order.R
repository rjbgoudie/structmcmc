# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Log likelihood of an order
#'
#' Evaluates the log likelihood of an order, using the efficient decompostion
#' used by Freidman and Koller (2003) (equation 8).
#'
#' @param order A vector length \code{numberOfNodes}, giving a permuation
#'   of \code{1:numberOfNodes}.
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
#' @references
#'   Friedman, N., & Koller, D. (2003). \emph{Being Bayesian about Network
#'   Structure. A Bayesian Approach to Structure Discovery in Bayesian
#'   Networks}. Machine Learning, 50, 95-125.
#'   \url{http://dx.doi.org/10.1023/A:1020249912095}
logOrderLikelihood <- function(order,
                               numberOfNodes,
                               nodesSeq,
                               scoresParents,
                               parentsTables,
                               allRows,
                               rowsThatContain){
  score <- vector("numeric", length = numberOfNodes)
  for (i in seq_along(order)){
    predecessors <- order[seq_len(i)]

    rows <- whichParentSetRows(node            = i,
                               nonDescendants  = predecessors,
                               numberOfNodes   = numberOfNodes,
                               allRows         = allRows,
                               rowsThatContain = rowsThatContain)
    score[i] <- logsumexp(scoresParents[[i]][rows])
  }
  sum(score)
}

#' Draw order proposal using flip-operator
#' 
#' Draws a proposal order, given a current order, using the single
#' flip-operator. See Grzegorcyck & Husmeier (2008) eq 12.
#' 
#' @param order A permutation of \code{nodesSeq}
#' @param nodesSeq The vector 1:nNodes(currentNetwork). (Supplied as an 
#'   argument for possible speed gain)
#' @return A new order
#' @export
orderFlipOperator <- function(order, nodesSeq){
  flip <- sample(nodesSeq, size = 2, replace = F)
  temp <- order[flip[1]]
  order[flip[1]] <- order[flip[2]]
  order[flip[2]] <- temp
  order
}


#' Order MCMC sampler for Bayesian Networks.
#'
#' Create a MCMC sampler for Bayesian Networks. The sampler samples Bayesian
#' Networks (ie models).
#'
#' @param data The data.
#' @param initial An object of class 'bn'. The starting value of the
#'                       MCMC.
#' @param orderPrior A function that returns the prior of an order.
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
BNOrderSampler <- function(data,
                           initial            = 1:ncol(data),
                           prior              = priorUniform(),
                           orderPrior         = orderPriorUniform,
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
  stopifnot(ncol(as.matrix(data)) ==   length(initial),
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

  initialDummy <- empty(numberOfNodes, "bn")
  constraint <- setupConstraint(constraint, initial = initialDummy)
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

  currentOrder <- initial
  rowsThatContain <- getRowsThatContain(numberOfNodes,
                                        parentsTables,
                                        maxNumberParents)

  allRows <- lapply(nodesSeq, function(node){
    seq_len(nrow(parentsTables[[node]]))
  })

  # Set up internal counters and logs etc
  nSteps <- 0

  sampler <- function(x,
                      verbose = F,
                      returnDiagnostics = F,
                      debugAcceptance = F,
                      returnTape = F,
                      burnin = 0){

    nBurnin <<- burnin
    nSteps <<- nSteps + 1

    # order MCMC
    proposalOrder <- orderFlipOperator(currentOrder, nodesSeq)

    logLikNewOrder <- logOrderLikelihood(order           = proposalOrder,
                                         numberOfNodes   = numberOfNodes,
                                         nodesSeq        = nodesSeq,
                                         scoresParents   = scoresParents,
                                         parentsTables   = parentsTables,
                                         allRows         = allRows,
                                         rowsThatContain = rowsThatContain)

    logLikOldOrder <- logOrderLikelihood(order           = currentOrder,
                                         numberOfNodes   = numberOfNodes,
                                         nodesSeq        = nodesSeq,
                                         scoresParents   = scoresParents,
                                         parentsTables   = parentsTables,
                                         allRows         = allRows,
                                         rowsThatContain = rowsThatContain)

    logAccProb <- logLikNewOrder -
                  logLikOldOrder +
                  log(orderPrior(proposalOrder)) -
                  log(orderPrior(currentOrder))

    logp <- log(runif(1, min = 0, max = 1))

    if (logAccProb >= 0 || logp < logAccProb){
      # if ACCEPTING the proposal
      currentOrder <<- proposalOrder
    } else {
      # if REJECTING the proposal
    }

    # finish off

    # if (isTRUE(keepTape)) updateTape(nSteps, currentNetwork)
    if (isTRUE(debugAcceptance)) browser()

    # updateET(currentNetwork, nSteps, nBurnin)
    # updateStatistics(currentNetwork, nSteps, nBurnin)

    currentOrder
  }
  class(sampler) <- c("sampler", "function")
  sampler
}

#' Sample a DAG given an order (weighted)
#' 
#' Sample a DAG given an order, using the parent weights. See Grzegorcyck & 
#' Husmeier (2008), eq 15.
#' 
#' @param order A vector length \code{numberOfNodes}, giving a permuation
#'   of \code{1:numberOfNodes}.
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
dagGivenOrder <- function(order,
                          numberOfNodes,
                          nodesSeq,
                          scoresParents,
                          parentsTables,
                          allRows,
                          rowsThatContain){
  out <- empty(numberOfNodes, "bn")
  scoreGivenOrder <- vector("numeric", numberOfNodes)
  for (i in seq_along(order)){
    predecessors <- order[seq_len(i)]
    node <- order[i]
    rows <- whichParentSetRows(node            = node,
                               nonDescendants  = predecessors,
                               numberOfNodes   = numberOfNodes,
                               allRows         = allRows,
                               rowsThatContain = rowsThatContain)
    scores <- scoresParents[[node]][rows]
    scoresNormalised <- exp(scores - logsumexp(scores))

    # sample a new parent set, according to the condtional probability
    samp <- sample.int(length(scores),
                       size = 1,
                       prob = scoresNormalised)
    scoreGivenOrder[i] <- scoresNormalised[samp]
    new <- parentsTables[[node]][rows[samp], ]
    out[[node]] <- new[!is.na(new)]
  }
  attr(out, "scoreGivenOrder") <- prod(scoreGivenOrder)
  out
}

orderPriorUniform <- function(x){
  1
}

orderPriorUniformStructures <- function(x){
  1/numberDAGsGivenOrder(x)
}

#' Is a BN consistent with an order?
#'
#' Tests whether a graph is consistent with a node ordering.
#'
#' @param x A \code{bn}.
#' @param order A vector length \code{numberOfNodes}, giving a permuation
#'   of \code{1:numberOfNodes}.
#' @return A logical indicating whether the BN is consistent with the
#'   supplied ordering.
#' @export
isConsistentWithOrder <- function(x, order){
  stopifnot(inherits(x, "bn"))
  is.ok <- T
  for (node in seq_along(order)){
    predecessors <- order[seq_len(which(order == node))]
    if (any(!x[[node]] %in% predecessors)){
      is.ok <- F
    }
  }
  is.ok
}

#' Number of BNs consistent with an order.
#'
#' Computes the number of BNs from a list of BNs that are consistent with an
#' order
#'
#' @param fam A \code{parental.list}
#' @param order A vector length \code{numberOfNodes}, giving a permuation
#'   of \code{1:numberOfNodes}.
#' @return A numeric. The number of BNs consistent with the order
#' @export
numberDAGsGivenOrder <- function(fam, order){
  length(fam[sapply(fam, isConsistentWithOrder, order = order)])
}

numberOrdersGivenDAG <- function(x){
  numberOfNodes <- length(x)
  nodesSeq <- seq_along(x)
  library(gtools)
  p <- permutations(3, 3, nodesSeq)

  ok <- vector("logical", nrow(p))
  for (i in 1:nrow(p)){
    order <- p[i, ]
    ok[i] <- isConsistentWithOrder(x, order)
  }
  sum(ok)
}

ellisWong <- function(orders, sampler, epsilon = 0.05){
  numberOfNodes <- get("numberOfNodes", env = environment(sampler))
  nodesSeq <- get("nodesSeq", env = environment(sampler))
  scoresParents <- get("scoresParents", env = environment(sampler))
  parentsTables <- get("parentsTables", env = environment(sampler))
  allRows <- get("allRows", env = environment(sampler))
  rowsThatContain <- get("rowsThatContain", env = environment(sampler))

  samples <- list()

  for (order in orders){
    cat(order, "\n")
    thisSamples <- parental.list()
    thisScores <- c()
    uniqueScores <- c()
    uniqueSamples <- parental.list()
    enough <- F
    i <- 1

    while (!enough){
      new <- dagGivenOrder(order,
                          numberOfNodes,
                          nodesSeq,
                          scoresParents,
                          parentsTables,
                          allRows,
                          rowsThatContain)
      thisSamples <- c(thisSamples, parental.list(new))
      class(thisSamples) <- "parental.list"
      thisScores <- c(thisScores, attr(new, "scoreGivenOrder"))
      if (!is.element2(new, uniqueSamples)){
        uniqueSamples <- c(uniqueSamples, parental.list(new))
        uniqueScores <- c(uniqueScores, attr(new, "scoreGivenOrder"))
      }
      i <- i + 1
      if (i > 1000){
        enough <- sampledEnough(uniqueScores, scores, epsilon)
      }
      if (i > 10000){
        warning("not finishing")
      }
    }
    samples <- c(samples, thisSamples)
  }
}

sampledEnough <- function(uniqueScores, scores, epsilon){
  sorted <- sort(uniqueScores, decreasing = T)
  uniqueScoresSeq <- seq_len(length(uniqueScores) - 1)
  ratios <- sorted[uniqueScoresSeq + 1]/sorted[uniqueScoresSeq]
  alpha <- max(ratios, na.rm = T)

  k <- log(epsilon * (1 - alpha) / max(uniqueScores))/log(alpha)
  k <- ceiling(k)

  numerator <- log(epsilon) - log(k)

  if (length(uniqueScores) >= k){
    denominator <- log(1 - uniqueScores[k])

    n <- length(scores)
    if (n >= numerator/denominator){
      T
    } else {
      F
    }
  } else {
    F
  }
}
