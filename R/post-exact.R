# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' BN Exact Posterior.
#'
#' method description
#'
#' @param bnspace ...
#' @param logScore ...
#' @param data ...
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
#' @export
#' @seealso \code{\link{bnpostmcmc}}, \code{\link{top.bnpost}}, 
#'   \code{\link{map.bnpost}}, \code{\link{gp.bnpost}},
#'   \code{\link{entropy.bnpost}}, \code{\link{tp.bnpost}},
#'   \code{\link{eht.bnpost}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#' 
#' bnspace <- enumerateBNSpace(2)
#' lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
#' post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)
#' 
#' gp(post)
#' ep(post)
bnpost <- function(bnspace, logScore, data, logScoreFUN = logScoreMultDir){
  stopifnot(class(bnspace)     == c("bn.list", "parental.list"),
            class(logScore)    == "numeric",
            class(logScoreFUN) == "function")

  out <- new.env()
  out[["type"]] <- "exact"
  out[["bnspace"]] <- bnspace
  out[["data"]] <- data
  out[["logScore"]] <- logScore
  out[["logScoreFUN"]] <- logScoreFUN
  class(out) <- "bnpost"
  out
}

#' Top graph from BN Posterior.
#'
#' method description
#'
#' @param x ...
#' @param head ...
#' @param ... Further arguments (unused)
#' @S3method top bnpost
#' @method top bnpost
#' @seealso \code{\link{top}}, \code{\link{map.bnpost}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#' 
#' bnspace <- enumerateBNSpace(2)
#' lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
#' post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)
#' 
#' top(post)
top.bnpost <- function(x, head = 10, ...){
  logScore <- x$logScore
  names(logScore) <- as.character(x$bnspace)
  rk <- rank(-logScore, ties.method = "min") <= head
  top <- sort(logScore[rk], dec = T)
  ids <- names(top)

  out <- as.parental(ids)
  out <- lapply(out, function(bn){
    class(bn) <- c("bn", "parental")
    bn
  })
  if (length(ids) > 1){
    class(out) <- c("bn.list", "parental.list")
  }
  out
}

#' Maximum aposteriori graph.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments (unused)
#' @S3method map bnpost
#' @method map bnpost
#' @seealso \code{\link{map}}, \code{\link{top.bnpost}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#' 
#' bnspace <- enumerateBNSpace(2)
#' lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
#' post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)
#' 
#' map(post)
map.bnpost <- function(x, ...){
  top(x, head = 1)
}

#' Posterior graph probaiblities.
#'
#' method description
#'
#' @param x ...
#' @param logNetworkPriors ...
#' @param log ...
#' @param pretty ...
#' @param ... Further arguments (unused)
#' @S3method gp bnpost
#' @method gp bnpost
#' @seealso \code{\link{gp}}, \code{\link{ep.bnpost}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#' 
#' bnspace <- enumerateBNSpace(2)
#' lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
#' post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)
#' 
#' gp(post)
#' ep(post)
gp.bnpost <- function(x, logNetworkPriors, log = F, pretty = F, ...){
  logScore <- x$logScore
  family <- x$bnspace

  if (missing(logNetworkPriors)){
    logNetworkPriors <- log(rep(1/length(family), length(family)))
  }

  normalisingConstant <- logsumexp(logScore + logNetworkPriors)
  out <- logScore + logNetworkPriors - normalisingConstant

  nm <- as.character(family, pretty)
  names(out) <- nm

  class(out) <- "gp"

  if (log){
    out
  }
  else {
    exp(out)
  }
}

#' Posterior edge probabilities.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments (unused)
#' @S3method ep bnpost
#' @method ep bnpost
#' @seealso \code{\link{ep}}, \code{\link{gp.bnpost}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#' 
#' bnspace <- enumerateBNSpace(2)
#' lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
#' post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)
#' 
#' gp(post)
#' ep(post)
ep.bnpost <- function(x, ...){
  logScore <- x$logScore
  family <- x$bnspace

  numberOfNodes <- length(family[[1]])
  nodesSeq <- seq_len(numberOfNodes)
  true <- matrix(0, numberOfNodes, numberOfNodes)
  logNetworkPriors <- log(rep(1/length(family), length(family)))
  normalisingConstant <- logsumexp(logScore + logNetworkPriors)

  for (head in nodesSeq){
    unlistedfamily <- unlist(family, rec = F)
    graphsFromHead <- unlistedfamily[seq(from = head,
                                         to = length(family) * numberOfNodes,
                                         by = numberOfNodes)]

    for (tail in nodesSeq){
      # get the graphs which have an edge tail --> head
      whichGraphs <- which(unlist(lapply(graphsFromHead, function(parents){
        if (tail %in% parents){
          T
        }
        else {
          F
        }})))

      if (length(whichGraphs) > 0){
        # sum up their scores
        true[tail, head] <- exp(logsumexp(logScore[whichGraphs] +
                                logNetworkPriors[whichGraphs]) -
                                normalisingConstant)
      }
    }
  }
  out <- true
  class(out) <- c("ep", "matrix")
  out
}

#' Entropy.
#'
#' method description
#'
#' @param x ...
#' @param logNetworkPriors ...
#' @param ... Further arguments (unused)
#' @S3method entropy bnpost
#' @method entropy bnpost
#' @seealso \code{\link{entropy}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#' 
#' bnspace <- enumerateBNSpace(2)
#' lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
#' post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)
#' 
#' entropy(post)
entropy.bnpost <- function(x, logNetworkPriors, ...){
  logScore <- x$logScore
  family <- x$bnspace
  if (missing(logNetworkPriors)){
    logNetworkPriors <- log(rep(1/length(family), length(family)))
  }

  #lgp <- graphProbs(logScore, family, logNetworkPriors, log = T)
  lgp <- gp(x, logNetworkPriors, log = T)

  ### dubious due to numerical overflow
  warning("worry about numerical overflow")
  -sum(exp(lgp) * lgp)
}

#' Transition probabilities.
#' 
#' Computes the matrix transition probabilities for the specified sampler.
#'
#' @param x An object of class "bnpost"
#' @param sampler Which sampler to use. Only "mh" for Metropolis-Hastings
#'               is implemented
#' @param allowFlips A logical of length 1, specifying whether the sampler is
#'               allowed to reverse the direction of single edges?
#' @param verbose ...
#' @param ... Further arguments (unused)
#' @return A matrix of transition probabilities.
#' @S3method tp bnpost
#' @method tp bnpost
#' @seealso \code{\link{tp}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#' 
#' bnspace <- enumerateBNSpace(2)
#' lsmd <- logScoreMultDir(bnspace, data = dat, hyperparameters = "qi")
#' post <- bnpost(bnspace = bnspace, logScore = lsmd, data = dat)
#' 
#' tp(post)
tp.bnpost <- function(x, sampler = "mh", allowFlips = T, verbose = F, ...){
  stopifnot(class(x) == "bnpost",
            sampler == "mh",
            class(allowFlips) == "logical",
            length(allowFlips) == 1)

  bnsp <- x$bnspace
  ng <- length(bnsp)
  p <- ncol(x$data)
  lsmd <- x$logScore

  nm <- matrix(NA, ng, ng)
  ngSeq <- seq.int(ng)
  if (verbose){
    cat("Computing SHD distances:")
  }
  for (i in ngSeq){
    if (verbose) cat("\n", i, ":")
    for (j in ngSeq){
      if (verbose) cat(j, ", ")
      nm[i, j] <- numberOfMovesBetweenIgnoringCycles(bnsp[[i]],
                                                     bnsp[[j]],
                                                     allowFlips = allowFlips)
    }
  }
  if (verbose) cat("\n Finishing off")
  nm0 <- ifelse(nm == 1, 1, 0)
  ns <- rowSums(nm0)

  # do this on the log scale
  alpha <- outer(-lsmd, lsmd, FUN = "+") + outer(log(ns), -log(ns), FUN = "+")
  alpha <- ifelse(alpha > 0, 1, exp(alpha))
  alpha <- outer(1/ns, rep(1, ng)) * alpha
  alpha <- alpha * nm0
  diag(alpha) <- 1 - rowSums(alpha)
  alpha
}

#' Hitting probability.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{eht}}, \code{\link{conductance}}, \code{\link{tp}}
hp <- function(x, ...){
  UseMethod("hp")
}

#' Expected hitting time.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{eht.bnpost}}
eht <- function(x, ...){
  UseMethod("eht")
}

#' Expected hitting time.
#' 
#' Computes the expected hitting times to the posterior modal graph from
#' the top 'head' other graphs.
#'
#' See Norris (1998) Markov Chains. Cambridge. Theorem 1.3.5, p17.
#'
#' @param x An object of class "bnpost".
#' @param head The expected hitting time from the top 'head' ranked graphs,
#'         with respect to the posterior will be returned. All are computed,
#'         so this setting has no effect on the run-time.
#' @param tp Optionally provide the transition probability matrix.
#' @param ... Further arguments (unused)
#'
#' @return A vector of the expected hitting times.
#' @S3method eht bnpost
#' @method eht bnpost
#' @seealso \code{\link{eht}}, \code{\link{tp}}, \code{\link{hp}}
eht.bnpost <- function(x, head = 5, tp = NULL, ...){
  stopifnot(class(x) == "bnpost",
            class(head) %in% c("numeric", "integer"),
            length(head) == 1)

  if (is.null(tp)){
    alpha <- tp(x)
  } else {
    alpha <- tp
  }
  diag(alpha) <- (1 - rowSums(alpha)) - 1

  ng <- nrow(alpha)
  o <- order(x$logScore, decreasing = T)
  whichMode <- o[1]
  alpha[whichMode, ] <- 0
  alpha[whichMode, whichMode] <- 1
  vec <- rep(-1, ng)
  vec[whichMode] <- 0

  rk <- rank(-x$logScore, ties.method = "min")

  if (is.finite(head)){
    wh <- c()
    for (i in seq_len(head)){
      wh <- c(wh, which(rk == i))
    }
  }
  else {
    wh <- T
  }
  solve(alpha, vec)[wh]
}

#' Transition probabilities.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{eht}}, \code{\link{hp}}, \code{\link{conductance}}
tp <- function(x, ...){
  UseMethod("tp")
}

#' Conductance.
#'
#' method description
#'
#' @param x ...
#' @param states ...
#' @export
#' @seealso \code{\link{eht}}, \code{\link{hp}}, \code{\link{tp}}
conductance <- function(x, states){
  # not complete nor correct!
  stop("not complete nor correct")
  stopifnot(class(x) == "bnpost",
            "bn.list" %in% class(states) | "bn" %in% class(states))
  scores <- gp(x)
  Pxy <- tp(x)
  S <- which2(states, x$bnspace)
  bnSeq <- seq_along(x$bnspace)

  out <- 0
  for (i in S){
    for (j in setdiff(bnSeq, S)){
      out <- out + scores[i] * Pxy[i, j]
    }
  }
  scoreS <- sum(scores[S])

  out/min(scoreS, 1 - scoreS)
}
