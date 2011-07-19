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

#' Local Normal-inverse-gamma (with g-prior) Log marginal likelihood.
#' 
#' Compute the LOCAL log marginal likelihood of the supplied
#' Bayesian Networks. ie the contribution to the log marginal liklihood from
#' one individual node.
#' 
#' Let \eqn{X}{X} be a data matrix with a number of predictors (in columns),
#' and \eqn{y}{y} be an response variable, and that \eqn{n}{n} observations
#' are available for each. For a graph \eqn{G}{G} (since this is local score
#' this is equivalent to an indicator vector), the model used is takes the
#' form
#' \eqn{y = \phi_{G} \beta + \epsilon}{y = phi_G * beta + epsilon}
#' with \eqn{\epsilon \sim N(0, \sigma^{2} I)}{epsilon ~ N(0, sigma^{2} I)}.
#' Note that the data needs to be standardised (zero-meaned).
#' 
#' The design matrix \eqn{\phi_{G}}{phi_{G}} is a column of 1s, and then
#' columns corresponding to each of the parents of the node. No cross-terms
#' are included.
#'
#' The prior used factorises as
#' \eqn{p(\beta, \sigma) = p(\beta \mid \sigma)p(\sigma)}{
#'     p(beta, sigma) = p(beta | sigma)p(sigma)}, 
#' The variance has an uninformative, scale invariant Jeffrey's prior 
#' \eqn{p(\sigma) = 1/\sigma^{2}}{p(sigma) = 1/sigma^2}, and the coefficients
#' have a zero-mean Normal prior (a Zellner g-prior), with \eqn{g = n}, so
#' that
#' \eqn{\beta \mid \sigma
#'      \sim
#'      N(0, g \sigma^{2}\left(\phi_{G}^{\prime}\phi_{G}\right)^{-1})}{
#'      beta | sigma ~ N(0, g * sigma^2 * (phi'_G phi_G)^-1)}
#' 
#' The above specification gives the following marginal likelihood.
#' 
#' \deqn{p(y \mid G)
#'       \propto
#'       (1 + n)^{-(\eta + 1)/2}
#'       \left(X^{T} X - \frac{n}{n + 1} X^{T}
#'        \phi_{G}(\phi^{T}_{G}\phi_{G})^{-1}\phi^{T}_{G}
#'        X\right)^{-\frac{n}{2}}
#'      }{
#'       P(y | G)
#'       propto
#'       (1 + n)^(-(eta + 1)/2) *
#'       (X' * X - (n/(n + 1)) * X' * 
#'        phi_G * (phi'_G * phi_G)^(-1) * phi_G *
#'        X)^(-n/2)
#'      }
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
#' @references
#'   Nott, D. J., & Green, P. J. (2004). \emph{Bayesian Variable Selection
#'   and the Swendsen-Wang Algorithm}. Journal of Computational and Graphical
#'   Statistics, 13, 141-157.
#'   \url{http://dx.doi.org/10.1198/1061860042958}
#'
#'   Smith, M., & Kohn, R. (1996). \emph{Nonparametric Regression using
#'   Bayesian Variable Selection}. Journal of Econometrics, 75, 317-343.
#'   \url{http://dx.doi.org/10.1016/0304-4076(95)01763-1}.
#' 
#'   Geiger, D., & Heckerman, D. (1994). \emph{Learning Gaussian Networks}.
#'   Proceedings of the 10th Conference Annual Conference on Uncertainty in
#'   Artificial Intelligence (UAI-94), 235-240.
#'   \url{http://uai.sis.pitt.edu/displayArticleDetails.jsp?mmnu=1&smnu=2&
#'   article_id=509&proceeding_id=10}
#' @export
#' @seealso \code{\link{logScoreNormal}},
#'   \code{\link{logScoreNormalOffline}},
#'   \code{\link{logScoreNormalIncremental}}.
localLogScoreNormal <- function(node,
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
    # phiQR <- qr(phi)
    # mx <- crossprod(Xi) - n/(n + 1) * crossprod(crossprod(qr.Q(phiQR), Xi))
    out <- -(eta + 1)/2 * log(1 + n) - n/2 * log(mx)
    (cache[[id]] <- out)
  }
}

#' Normal-inverse-gamma (with g-prior) Log marginal likelihood.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{logScoreNormal.bn}}
#' @examples
#' data <- cbind(c(-10, -2), c(-11, -4))
#' net <- bn(integer(0), 1)
#' logScoreNormal(net, data)
logScoreNormal <- function(x, ...){
  UseMethod("logScoreNormal")
}

#' Normal-inverse-gamma (with g-prior) Log marginal likelihood.
#' 
#' Compute the log marginal likelihood of the supplied Bayesian Network.
#'
#' The data is scored as continuous, using a form of the Normal Prior.
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
#' @S3method logScoreNormal bn
#' @method logScoreNormal bn
#' @seealso \code{\link{logScoreNormal}},
#'   \code{\link{logScoreNormal.bn.list}}
#' @examples
#' data <- cbind(c(-10, -2), c(-11, -4))
#' net <- bn(integer(0), 1)
#' logScoreNormal(net, data)
logScoreNormal.bn <- function(x,
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
  logScoreParameters <- logScoreNormalPrepare(data               = data,
                                              logScoreParameters = list(),
                                              checkInput = F)
  sum(.Internal(unlist(lapply(nodeSeq, function(i){
    localLogScoreNormal(node               = i,
                        parents            = x[[i]],
                        logScoreParameters = logScoreParameters,
                        cache              = cache)
  }), F, F)))
}

#' Normal-inverse-gamma (with g-prior) Log marginal likelihood (offline).
#' 
#' Compute the log marginal likelihood of the supplied Bayesian Network.
#'
#' This function is an alternative interface to logScoreNormal.
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
#' @seealso \code{\link{logScoreNormal}},
#'   \code{\link{logScoreNormalIncremental}}
logScoreNormalOffline <- function(x,
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
    localLogScoreNormal(node               = head,
                        parents            = x[[head]],
                        logScoreParameters = logScoreParameters,
                        cache              = cache,
                        checkInput         = F)
  }), F, F)))
}

#' Normal-inverse-gamma (with g-prior) Log marginal likelihood.
#' 
#' Compute the log marginal likelihood of the supplied Bayesian Networks.
#'
#' The data is scored as continuous, using a form of the Normal Prior.
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
#' @S3method logScoreNormal bn.list
#' @method logScoreNormal bn.list
#' @seealso \code{\link{logScoreNormal.bn}}, \code{\link{logScoreNormal}}
logScoreNormal.bn.list <- function(x,
                                   data,
                                   cache = new.env(hash = TRUE),
                                   ...){
  stopifnot(class(data)                     ==   "matrix",
            all(unlist(lapply(data, class)) %in% c("numeric", "integer")),
            "bn.list"                       %in% class(x))
  unlist(lapply(x, function(bn){
    logScoreNormal(bn, data, cache)
  }))
}

#' Internal functions.
#' 
#' Convert a data frame to the appropriate format for the fast/incremental
#' logScoreNormal functions, and return as part of the logScoreParameters
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
#' @seealso \code{\link{logScoreNormal}}
logScoreNormalPrepare <- function(data,
                                  logScoreParameters,
                                  checkInput = T){
  if (isTRUE(checkInput)){
    stopifnot(class(data)                     == "data.frame",
              all(unlist(lapply(data, class)) == "factor"))
  }
  nl <- apply(data, 2, nunique)
  modifyList(logScoreParameters, val = list(data = data, nl = nl))
}

#' Normal-inverse-gamma (with g-prior) Log marginal likelihood (online).
#' 
#' Compute the difference in log marginal likelihood of the supplied
#' Bayesian Networks, quickly.
#'
#' This is a fast, incremental version of logScoreNormal.
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
#' @seealso \code{\link{logScoreNormal}},
#'   \code{\link{logScoreNormalOffline}}
logScoreNormalIncremental <- function(currentBN,
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
    localLogScoreNormal(node               = head,
                        parents            = proposalBN[[head]],
                        logScoreParameters = logScoreParameters,
                        cache              = cache,
                        checkInput         = F) -
    cache[[fastid(c(head, currentBN[[head]]))]]
  }), F, F)))
}
