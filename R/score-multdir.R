# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Fast ID of a graph.
#'
#' method description
#'
#' @param id ...
#' @export
fastid <- function(id){
  .Internal(paste(list(id), "", ","))
}

#' Checks validity.
#'
#' Checks whether the supplied hyperparameters parameter is valid.
#'
#' @param x The hyperparameters to test.
#'
#' @return A logical of length 1.
#' @S3method is.valid hyp
#' @method is.valid hyp
#' @seealso \code{\link{logScoreMultDir}}
is.valid.hyp <- function(x){
  tryCatch({
    all(class(x)  ==   "character",
        length(x) ==   1,
        x         %in% c("qi", "one", "point9"))
    },
    error = function(e) F)
}

#' Local Multinomial-Dirichlet Log marginal likelihood.
#'
#' Compute the LOCAL log marginal likelihood of the supplied
#' Bayesian Networks. ie the contribution to the log marginal liklihood from
#' one individual node.
#'
#' The data must be discrete. The conditional distributions of each
#' random variable, conditional on its parents are assumed to be
#' multinomial, with Dirichlet priors for the parameters.
#'
#' The notation here roughly follows
#' Mukherjee and Speed (2008) Network inference using informative priors.
#' PNAS 105 (38) 14313-14318, doi: 10.1073/pnas.0802272105
#'
#' @param node A numeric vector of length 1. The node to compute the local
#'   log score for.
#' @param parents A numeric vector. The parents of node.
#' @param logScoreParameters A list with the following components:
#'   \describe{
#'     \item{data}{A matrix (NOT data.frame), with columns being integers
#'                 in the range 0, 1, 2, ....  giving the values of each
#'                 random variable.
#'                 **** The integers MUST start  numbering at 0 NOT 1    ****}
#'     \item{nl}{A numeric vector of length ncol(data), specifying the number
#'               of levels that each random variable takes.}
#'     \item{hyperparameters}{A character vector of length one.
#'                            Either "qi", "one", or "point9"}
#'   }
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param checkInput A logical of length 1, specifying whether to check the
#'   inputs to the function.
#' @return A numeric vector of length 1, giving the log marginal likelihood.
#'   The environment 'cache' will also be updated because its scope is
#'   global.
#' @export
#' @seealso \code{\link{logScoreMultDir}},
#'   \code{\link{logScoreMultDirIncremental}},
#'   \code{\link{logScoreMultDirOffline}}
localLogScoreMultDir <- function(node,
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
              is.valid.hyp(logScoreParameters$hyperparameters) ||
                is.null(logScoreParameters$hyperparameters),
              class(cache)                   ==   "environment")
  }
  if (is.null(logScoreParameters$hyperparameters)){
    logScoreParameters$hyperparameters <- "qi"
  }
  whichCols <- c(node, parents)
  id <- fastid(whichCols)
  cacheRecord <- cache[[id]]
  if (!is.null(cacheRecord)){
    cacheRecord
  }
  else {
    d <- logScoreParameters$nl[whichCols]
    if (length(whichCols) > 1){
      N_prime <- switch(logScoreParameters$hyperparameters,
        qi = 1/prod(d[-1]),
        one = 1,
        point9 = 0.9,
        1/prod(d[-1]) # default is qi
      )
    }
    else {
      N_prime <- switch(logScoreParameters$hyperparameters,
          qi = 1,
          one = 1,
          point9 = 0.9
        )
    }
    N_prime_marginals <- d[1] * N_prime
    pd <- cumprod(c(1, d))
    sel <- length(d) + 1
    if (pd[sel] > 2147483647L){
      stop("Cannot compute this, due to too many parents")
    }
    # tcrossprod would be faster, but there are bugs in < R 2.9.1
    bin <- tcrossprod(pd[-sel], logScoreParameters$data[, whichCols]) + 1
    # bin <- pd[-sel] %*% t(data[, whichCols]) + 1
    N <- .C("R_tabulate",
            as.integer(bin),
            as.integer(length(bin)),
            as.integer(pd[sel]),
            ans = integer(pd[sel]),
            PACKAGE = "base")$ans
    (cache[[id]] <- sum(lgamma(N_prime + N) -
                        lgamma(N_prime)) -
                    sum(lgamma(N_prime_marginals +
                      colSums(matrix(N,
                              nrow = logScoreParameters$nl[whichCols[1]]))) -
                        lgamma(N_prime_marginals)))
  }
}

#' Multinomial-Dirichlet Log marginal likelihood.
#'
#' method description
#'
#' @param x A \code{bn}.
#' @param ... Further arguments, passed to method
#' @export
#' @seealso \code{\link{logScoreMultDir.bn}},
#'   \code{\link{logScoreMultDir.bn.list}},
#'   \code{\link{logScoreMultDirOffline}},
#'   \code{\link{logScoreMultDirIncremental}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
#' x2 <- factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
#' x3 <- factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
#' data <- data.frame(x1 = x1, x2 = x2,  x3 = x3)
#' logScoreMultDir(bn(integer(0), 1, 2), data)
logScoreMultDir <- function(x, ...){
  UseMethod("logScoreMultDir")
}

#' Compute the log marginal likelihood of the supplied Bayesian Network.
#'
#' The data must be discrete. The conditional distributions of each
#' random variable, conditional on its parents are assumed to be
#' multinomial, with Dirichlet priors for the parameters.
#'
#' The notation here roughly follows
#' Mukherjee and Speed (2008) Network inference using informative priors.
#' PNAS 105 (38) 14313-14318, doi: 10.1073/pnas.0802272105
#'
#' @param x An object of class "bn". The Bayesian Network by for which the
#'   marginal likelihood is computed.
#' @param data A data.frame, with columns being factors giving the values of
#'   each random variable.
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param hyperparameters A character vector of length one. Either "qi",
#'   "one", or "point9"
#' @param checkInput A logical of length 1, specifying whether to check the
#'   inputs to the function.
#' @param ... Further arguments, currently unused
#' @return A numeric vector of length 1, giving the log marginal likelihood.
#'   The environment 'cache' will also be updated because its scope is
#'   global.
#' @S3method logScoreMultDir bn
#' @method logScoreMultDir bn
#' @seealso \code{\link{logScoreMultDir}},
#'   \code{\link{logScoreMultDir.bn.list}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
#' x2 <- factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
#' x3 <- factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
#' data <- data.frame(x1 = x1, x2 = x2,  x3 = x3)
#' logScoreMultDir(bn(c(), c(1), c(2)), data)
logScoreMultDir.bn <- function(x,
                               data,
                               cache           = new.env(hash = T),
                               hyperparameters = "qi",
                               checkInput      = T,
                               ...){
  if (isTRUE(checkInput)){
    stopifnot("bn"                            %in% class(x),
              is.valid(x),
              class(data)                     ==   "data.frame",
              all(unlist(lapply(data, class)) ==   "factor"),
              class(cache)                    ==   "environment",
              is.valid.hyp(hyperparameters))
  }
  p <- nNodes(x)
  nodeSeq <- seq_len(p)
  logScoreParameters <- logScoreMultDirPrepare(data               = data,
                                               logScoreParameters =
                                                 list(hyperparameters =
                                                   hyperparameters),
                                               checkInput = F)
  sum(.Internal(unlist(lapply(nodeSeq, function(i){
    localLogScoreMultDir(node               = i,
                         parents            = x[[i]],
                         logScoreParameters = logScoreParameters,
                         cache              = cache)
  }), F, F)))
}

#' Multinomial-Dirichlet Log marginal likelihood (offline).
#'
#' Compute the log marginal likelihood of the supplied Bayesian Network.
#'
#' This function is an alternative interface to logScoreMultDir.
#' This interface is required by the MCMC sampler.
#'
#' @param x An object of class "bn". The Bayesian Network by for which the
#'   marginal likelihood is computed.
#' @param logScoreParameters A list with the following components:
#'   \describe{
#'     \item{data}{A matrix (NOT data.frame), with columns being integers
#'                 in the range 0, 1, 2, ....  giving the values of each
#'                 random variable.
#'                 **** The integers MUST start  numbering at 0 NOT 1    ****}
#'     \item{nl}{A numeric vector of length ncol(data), specifying the number
#'               of levels that each random variable takes.}
#'     \item{hyperparameters}{A character vector of length one.
#'                            Either "qi", "one", or "point9"}
#'   }
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param checkInput A logical of length 1, specifying whether to check the
#'   inputs to the function.
#' @return A numeric vector of length 1, giving the log marginal likelihood.
#'   The environment 'cache' will also be updated because its scope is
#'   global.
#' @export
#' @seealso \code{\link{logScoreMultDir}},
#'   \code{\link{logScoreMultDirIncremental}}
logScoreMultDirOffline <- function(x,
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
              is.valid.hyp(logScoreParameters$hyperparameters) ||
                is.null(logScoreParameters$hyperparameters),
              class(cache)                   ==   "environment")
  }
  sum(.Internal(unlist(lapply(seq_len(nNodes(x)), function(head){
    localLogScoreMultDir(node               = head,
                         parents            = x[[head]],
                         logScoreParameters = logScoreParameters,
                         cache              = cache,
                         checkInput         = F)
  }), F, F)))
}

#' Multinomial-Dirichlet Log marginal likelihood.
#'
#' Compute the log marginal likelihood of the supplied Bayesian Networks.
#'
#' The data must be discrete. The conditional distributions of each
#' random variable, conditional on its parents are assumed to be
#' multinomial, with Dirichlet priors for the parameters.
#'
#' The notation here roughly follows
#' Mukherjee and Speed (2008) Network inference using informative priors.
#' PNAS 105 (38) 14313-14318, doi: 10.1073/pnas.0802272105
#'
#' @param x An object of class "bn.list", the Bayesian Networks for which
#'   the marginal likelihood are computed.
#' @param data A data.frame, with columns being factors giving the values of
#'   each random variable.
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param hyperparameters  A character vector of length one. Either "qi",
#'   "one", or "point9"
#' @param verbose A logical of length 1. If true, a progress bar will
#'   be shown.
#' @param ... Further arguments (unused)
#' @return A numeric vector of length 1, giving the log marginal likelihood.
#'   The environment 'cache' will also be updated because its scope is
#'   global.
#' @S3method logScoreMultDir bn.list
#' @method logScoreMultDir bn.list
#' @seealso \code{\link{logScoreMultDir}}, \code{\link{logScoreMultDir.bn}},
#'   \code{\link{logScoreMultDirOffline}},
#'   \code{\link{logScoreMultDirIncremental}}
logScoreMultDir.bn.list <- function(x,
                                    data,
                                    hyperparameters = "qi",
                                    cache           = new.env(hash = T),
                                    verbose = F,
                                    ...){
  stopifnot(class(data)                     ==   "data.frame",
            all(unlist(lapply(data, class)) ==   "factor"),
            "bn.list"                       %in% class(x))
  if (verbose){
    progress <- txtProgressBar(max = length(x), style = 3)
     setTxtProgressBar(progress, 0)
     prog <- 0
  }
  out <- unlist(lapply(x, function(bn){
    if (verbose){
      prog <<- prog + 1
      setTxtProgressBar(progress, prog)
    }
    logScoreMultDir(bn, data, cache, hyperparameters = hyperparameters)
  }))
  if (verbose){
    close(progress)
  }
  out
}

#' Internal functions.
#'
#' Convert a data frame to the appropriate format for the fast/incremental
#' logScoreMultDir functions, and return as part of the logScoreParameters
#' list.
#'
#' In particular, the data is converted to a matrix, and the factor levels
#' taken integer values from 0, 1, .... ie not on 1, 2, 3.
#'
#' @param data A data.frame, with columns being factors giving the values of
#'    each random variable.
#' @param logScoreParameters A list, optionally containing other parameters.
#'   In particular, it may contain 'hyperparameters'.
#'   \describe{
#'     \item{data}{A matrix (NOT data.frame), with columns being integers
#'                 in the range 0, 1, 2, ....  giving the values of each
#'                 random variable.
#'                 **** The integers MUST start  numbering at 0 NOT 1    ****}
#'     \item{nl}{A numeric vector of length ncol(data), specifying the number
#'               of levels that each random variable takes.}
#'     \item{hyperparameters}{A character vector of length one.
#'                            Either "qi", "one", or "point9"}
#'   }
#' @param checkInput A logical of length 1, specifying whether to check
#'   the inputs to the function.
#' @return A list with the contents of logScoreParameters, with the following
#'   components added or altered:
#'   \describe{
#'     \item{data}{A matrix (NOT data.frame), with columns being integers
#'                 in the range 0, 1, 2, ....  giving the values of each
#'                 random variable.}
#'     \item{nl}{A numeric vector of length ncol(data), specifying the number
#'               of levels that each random variable takes.}
#'   }
#' @export
#' @seealso \code{\link{logScoreMultDir}}
logScoreMultDirPrepare <- function(data, logScoreParameters, checkInput = T){
  if (isTRUE(checkInput)){
    stopifnot(class(data)                     == "data.frame",
              all(unlist(lapply(data, class)) == "factor"))
  }
  nl <- sapply(data, nlevels)
  data <- data.matrix(data) - 1
  modifyList(logScoreParameters, val = list(data = data, nl = nl))
}

#' Multinomial-Dirichlet Log marginal likelihood (online).
#'
#' Compute the difference in log marginal likelihood of the supplied
#' Bayesian Networks, quickly.
#'
#' This is a fast, incremental version of logScoreMultDir.
#'
#' *** NOTE: THIS REQUIRES A MATRIX RATHER THAN A DATA FRAME ***
#'
#' The data must be discrete. The conditional distributions of each
#' random variable, conditional on its parents are assumed to be
#' multinomial, with Dirichlet priors for the parameters.
#'
#' The notation here roughly follows
#' Mukherjee and Speed (2008) Network inference using informative priors.
#' PNAS 105 (38) 14313-14318, doi: 10.1073/pnas.0802272105
#'
#' @param currentBN An object of class "bn".
#' @param proposalBN An object of class "bn".
#' @param heads A numeric vector, specifying which nodes have different
#'   parents in currentBN and proposalBN.
#' @param logScoreParameters A list with the following components:
#'   \describe{
#'     \item{data}{A matrix (NOT data.frame), with columns being integers
#'                 in the range 0, 1, 2, ....  giving the values of each
#'                 random variable.
#'                 **** The integers MUST start  numbering at 0 NOT 1    ****}
#'     \item{nl}{A numeric vector of length ncol(data), specifying the number
#'               of levels that each random variable takes.}
#'     \item{hyperparameters}{A character vector of length one.
#'                            Either "qi", "one", or "point9"}
#'   }
#' @param cache Optionally, provide an environment with cached local scores
#'   for this data.
#' @param checkInput A logical of length 1, specifying whether to check the
#'   inputs to the function.
#' @return logscore(proposalBN) - logscore(currentBN)
#' @export
#' @seealso \code{\link{logScoreMultDir}},
#'   \code{\link{logScoreMultDirOffline}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1, 0, 0, 1, 0, 1, 0))
#' x2 <- factor(c(0, 1, 0, 1, 0, 1, 1, 0, 1, 0))
#' x3 <- factor(c(0, 1, 1, 1, 0, 1, 1, 0, 1, 0))
#' data <- data.frame(x1 = x1, x2 = x2,  x3 = x3)
#'
#' n1 <- bn(numeric(0), 1, 2)
#' n2 <- bn(numeric(0), 1, numeric(0))
#'
#' logScoreMultDir(n1, data)
#' logScoreMultDir(n2, data)
#'
#' data <- data.frame(lapply(data, as.factor))
#' data <- data.matrix(data) - 1
#' cc <- new.env(hash = TRUE, size = 10000L)
#' nl <- apply(data, 2, function(i) length(unique(i)))
#' names(nl) <- seq_along(n1)
#' lsp <- list(data = data, nl = nl, hyperparameters = "qi")
#' logScoreMultDirIncremental(n1, n2, 3, lsp, cc)
logScoreMultDirIncremental <- function(currentBN,
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
              is.valid.hyp(logScoreParameters$hyperparameters) ||
                is.null(logScoreParameters$hyperparameters),
              class(cache)                   ==   "environment")
  }
  sum(.Internal(unlist(lapply(heads, function(head){
    localLogScoreMultDir(node               = head,
                         parents            = proposalBN[[head]],
                         logScoreParameters = logScoreParameters,
                         cache              = cache,
                         checkInput         = F) -
    cache[[fastid(c(head, currentBN[[head]]))]]
  }), F, F)))
}
