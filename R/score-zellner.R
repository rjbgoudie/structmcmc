# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Number of unique.
#'
#' method description
#'
#' @param x ...
#' @export
nunique <- function(x){
  length(unique(x))
}

#' Local Normal Log marginal likelihood.
#' 
#' Compute the LOCAL log marginal likelihood of the supplied
#' Bayesian Networks. ie the contribution to the log marginal liklihood from
#' one individual node.
#'
#' The data is scored as continuous, using a form of the Zellner Prior.
#'
#' @param node A numeric vector of length 1. The node to compute the local
#'   log score for.
#' @param parents A numeric vector. The parents of node.
#' @param logScoreParameters A list with the following components:
#'   \describe{
#'     \item{data}{A matrix with columns giving the values of each random
#'                 variable.}
#'     \item{nl}{A numeric vector of length nNodes(currentBN), specifying the
#'               number of levels that each random variable takes.}
#'   }
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param checkInput A logical of length 1, specifying whether to check the
#'   inputs to the function.
#' @return A numeric vector of length 1, giving the log marginal likelihood.
#'   The environment 'cache' will also be updated because its scope is
#'   global.
#' @export
#' @seealso \code{\link{logScoreZellner}},
#'   \code{\link{logScoreZellnerOffline}},
#'   \code{\link{logScoreZellnerIncremental}}.
localLogScoreZellner <- function(node,
                                 parents,
                                 logScoreParameters,
                                 cache,
                                 checkInput = T){
  if (isTRUE(checkInput)){
    stopifnot(class(node)                    %in% c("numeric", "integer"),
              length(node)                   ==   1,
              class(parents)                 %in% c("numeric", "integer"),
              is.list(logScoreParameters),
              class(logScoreParameters$data) ==   "matrix",
              class(logScoreParameters$nl)   %in% c("numeric", "integer"),
              class(cache)                   ==   "environment")
  }
  id <- fastid(c(node, parents))
  cacheRecord <- cache[[id]]
  if (!is.null(cacheRecord)){
    cacheRecord
  }
  else {
    eta <- length(parents)
    n <- nrow(logScoreParameters$data)
    Xi <- logScoreParameters$data[, node]
    phi <- matrix(data = c(rep(1, n), logScoreParameters$data[, parents]),
                  ncol = eta + 1,
                  nrow = n)
    # solve should throw an error on singularity, or ill-conditioning.
    mx <- crossprod(Xi) -
          n/(n + 1) * crossprod(Xi, phi) %*%
          solve(crossprod(phi)) %*% crossprod(phi, Xi)
    out <- -(eta + 1)/2 * log(1 + n) - n/2 * log(mx)
    (cache[[id]] <- out)
  }
}

#' Normal Log marginal likelihood.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{logScoreZellner.bn}}
logScoreZellner <- function(x, ...){
  UseMethod("logScoreZellner")
}

#' Normal Log marginal likelihood.
#' 
#' Compute the log marginal likelihood of the supplied Bayesian Network.
#'
#' The data is scored as continuous, using a form of the Zellner Prior.
#'
#' @param x An object of class "bn". The Bayesian Network by for which the
#'   marginal likelihood is computed.
#' @param data A matrix with columns giving the value of each random variable.
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param checkInput A logical of length 1, specifying whether to check the
#'   inputs to the function.
#' @param ... Further arguments, currently unused
#' @return A numeric vector of length 1, giving the log marginal likelihood.
#'   The environment 'cache' will also be updated because its scope is
#'   global.
#' @S3method logScoreZellner bn
#' @method logScoreZellner bn
#' @seealso \code{\link{logScoreZellner}},
#'   \code{\link{logScoreZellner.bn.list}}
logScoreZellner.bn <- function(x,
                               data,
                               cache      = new.env(hash = T),
                               checkInput = T,
                               ...){
  if (isTRUE(checkInput)){
    stopifnot("bn"                            %in% class(x),
              is.valid(x),
              class(data)                     ==   "matrix",
              all(unlist(lapply(data, class)) %in% c("numeric", "integer")),
              class(cache)                    ==   "environment")
  }
  p <- nNodes(x)
  nodeSeq <- seq_len(p)
  logScoreParameters <- logScoreZellnerPrepare(data               = data,
                                               logScoreParameters = list(),
                                               checkInput = F)
  sum(.Internal(unlist(lapply(nodeSeq, function(i){
    localLogScoreZellner(node               = i,
                         parents            = x[[i]],
                         logScoreParameters = logScoreParameters,
                         cache              = cache)
  }), F, F)))
}

#' Normal Log marginal likelihood (offline).
#' 
#' Compute the log marginal likelihood of the supplied Bayesian Network.
#'
#' This function is an alternative interface to logScoreZellner.
#' This interface is required by the MCMC sampler.
#'
#' @param x An object of class "bn". The Bayesian Network by for which the
#'   marginal likelihood is computed.
#' @param logScoreParameters A list with the following components:
#'   \describe{
#'     \item{data}{A matrix with columns giving the values of each random
#'                 variable.}
#'     \item{nl}{A numeric vector of length nNodes(currentBN), specifying the
#'               number of levels that each random variable takes.}
#'   }
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param checkInput A logical of length 1, specifying whether to check the
#'   inputs to the function.
#' @return A numeric vector of length 1, giving the log marginal likelihood.
#'   The environment 'cache' will also be updated because its scope is
#'   global.
#' @export
#' @seealso \code{\link{logScoreZellner}},
#'   \code{\link{logScoreZellnerIncremental}}
logScoreZellnerOffline <- function(x,
                                   logScoreParameters,
                                   cache      = new.env(hash = T),
                                   checkInput = T){
  if (isTRUE(checkInput)){
    stopifnot("bn"                           %in% class(x),
              is.valid(x),
              is.list(logScoreParameters),
              class(logScoreParameters$data) ==   "matrix",
              class(logScoreParameters$nl)   %in% c("numeric", "integer"),
              length(logScoreParameters$nl)  ==   nNodes(x),
              class(cache)                   ==   "environment")
  }
  sum(.Internal(unlist(lapply(seq_len(nNodes(x)), function(head){
    localLogScoreZellner(node               = head,
                         parents            = x[[head]],
                         logScoreParameters = logScoreParameters,
                         cache              = cache,
                         checkInput         = F)
  }), F, F)))
}

#' Normal Log marginal likelihood.
#' 
#' Compute the log marginal likelihood of the supplied Bayesian Networks.
#'
#' The data is scored as continuous, using a form of the Zellner Prior.
#'
#' @param x An object of class "bn.list", the Bayesian Networks for which
#'   the marginal likelihood are computed.
#' @param data A matrix, with columns giving the values of each random
#'   variable.
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param ... Further arguments, currently unused.
#'
#' @return A numeric vector of length 1, giving the log marginal likelihood.
#'   The environment 'cache' will also be updated because its scope is
#'   global.
#' @S3method logScoreZellner bn.list
#' @method logScoreZellner bn.list
#' @seealso \code{\link{logScoreZellner.bn}}, \code{\link{logScoreZellner}}
logScoreZellner.bn.list <- function(x,
                                    data,
                                    cache = new.env(hash = TRUE),
                                    ...){
  stopifnot(class(data)                     ==   "matrix",
            all(unlist(lapply(data, class)) %in% c("numeric", "integer")),
            "bn.list"                       %in% class(x))
  unlist(lapply(x, function(bn){
    logScoreZellner(bn, data, cache)
  }))
}

#' Internal functions.
#' 
#' Convert a data frame to the appropriate format for the fast/incremental
#' logScoreZellner functions, and return as part of the logScoreParameters
#' list.
#'
#' In particular, the data is converted to a matrix, and the factor levels
#' taken integer values from 0, 1, .... ie not on 1, 2, 3.
#'
#' @param data A data.frame, with columns being factors giving the values of
#'   each random variable.
#' @param logScoreParameters A list with the following components:
#'   \describe{
#'     \item{data}{A matrix with columns giving the values of each random
#'                 variable.}
#'     \item{nl}{A numeric vector of length nNodes(currentBN), specifying the
#'               number of levels that each random variable takes.}
#'   }
#' @param checkInput A logical of length 1, specifying whether to check the
#'   inputs to the function.
#' @return A list with the contents of logScoreParameters, with the following
#'   components added or altered:
#'   \describe{
#'     \item{data}{A matrix with columns giving the value of each random
#'                 variable.}
#'     \item{nl}{A numeric vector of length ncol(data), specifying the number
#'               of levels that each random variable takes.}
#'   }
#' @export
#' @seealso \code{\link{logScoreZellner}}
logScoreZellnerPrepare <- function(data,
                                   logScoreParameters,
                                   checkInput = T){
  if (isTRUE(checkInput)){
    stopifnot(class(data)                     == "data.frame",
              all(unlist(lapply(data, class)) == "factor"))
  }
  nl <- apply(data, 2, nunique)
  modifyList(logScoreParameters, val = list(data = data, nl = nl))
}

#' Normal Log marginal likelihood (online).
#' 
#' Compute the difference in log marginal likelihood of the supplied
#' Bayesian Networks, quickly.
#'
#' This is a fast, incremental version of logScoreZellner.
#'
#' The data is scored as continuous, using a form of the Zellner Prior.
#'
#' @param currentBN An object of class "bn".
#' @param proposalBN An object of class "bn".
#' @param heads A numeric vector, specifying which nodes have different
#'   parents in currentBN and proposalBN.
#' @param logScoreParameters A list with the following components:
#'   \describe{
#'     \item{data}{A matrix with columns giving the values of each random
#'                 variable.}
#'     \item{nl}{A numeric vector of length nNodes(currentBN), specifying the
#'               number of levels that each random variable takes.}
#'   }
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param checkInput A logical of length 1, specifying whether to check the
#'   inputs to the function.
#' @return logscore(proposalBN) - logscore(currentBN)
#' @export
#' @seealso \code{\link{logScoreZellner}},
#'   \code{\link{logScoreZellnerOffline}}
logScoreZellnerIncremental <- function(currentBN,
                                       proposalBN,
                                       heads,
                                       logScoreParameters,
                                       cache,
                                       checkInput = T){
  if (isTRUE(checkInput)){
    stopifnot("bn"                           %in% class(currentBN),
              is.valid(currentBN),
              "bn"                           %in% class(proposalBN),
              is.valid(proposalBN),
              nNodes(currentBN)              ==   nNodes(proposalBN),
              class(heads)                   %in% c("numeric", "integer"),
              is.list(logScoreParameters),
              class(logScoreParameters$data) ==   "matrix",
              class(logScoreParameters$nl)   %in% c("numeric", "integer"),
              length(logScoreParameters$nl)  ==   nNodes(currentBN),
              class(cache)                   ==   "environment")
  }
  sum(.Internal(unlist(lapply(heads, function(head){
    localLogScoreZellner(node               = head,
                         parents            = proposalBN[[head]],
                         logScoreParameters = logScoreParameters,
                         cache              = cache,
                         checkInput         = F) -
    cache[[fastid(c(head, currentBN[[head]]))]]
  }), F, F)))
}
