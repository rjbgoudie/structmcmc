# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Mukherjee Bioinformatics prior.
#'
#' Compute prior score P(M) for a 'bvsresponse'. The prior is flat from
#' 0 parents to k0 parents, then exponentially(lambda) decreaasing to kmax,
#' at which point it becomes zero.
#'
#' The prior originates in Mukherjee el al. (2009), and is equation (8) in
#' section 3.1.3.
#'
#' Mukherjee et al. Sparse combinatorial inference with an application in
#' cancer biology. Bioinformatics (2009) vol. 25 (2) pp. 265
#'
#' @param x A 'bvsresponse' object
#' @param k0 The level below which the prior stops being flat
#' @param kmax The level above which the prior is zero
#' @param lambda The exponential decay parameter
#' @return The prior score for the 'bvsresponse' x
#' @export
mukherjeeBioinformaticsPrior <- function(x, k0, kmax, lambda){
  stopifnot("bvsresponse" %in% class(x))

  numberOfParents <- length(x$parents)
  kmaxTest <- numberOfParents <= kmax
  if (kmaxTest){
    part <- min(0, k0 - numberOfParents)
    exp(lambda * part)
  }
  else {
    0
  }
}

#' A standard 'graph prior'.
#'
#' Returns a function that will evaluate the prior of a graph. The prior is
#'
#' f(x) = exp( - lambda * (edges(x) - edges(graph)))
#'
#' So the prior scores more highly those graphs \code{x} that include all the
#' edges in the graph \code{graph}; extra edges are penalised; the prior
#' is maximised (not uniquely) at the empty graph.
#'
#' @param graph The 'prior graph'. A \code{bn}.
#' @param lambda A weighting parameter. A numeric of length 1.
#' @return A function computes the prior score of the supplied graph. This
#'   This function is of a suitable form to be used as a prior.
#' @export
#' @seealso \code{\link{priorUniform}}
#' @examples
#' x1 <- factor(c("a", "a", "g", "c", "c", "a", "g", "a", "a"))
#' x2 <- factor(c(2, 2, 4, 3, 1, 4, 4, 4, 1))
#' x3 <- factor(c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
#' x <- data.frame(x1 = x1, x2 = x2, x3 = x3)
#'
#' priorgraph <- bn(c(), 1, 2)
#' prior <- priorGraph(priorgraph, 0.5)
#'
#' initial <- empty(3, "bn")
#' sampler <- BNSampler(data = x, initial = initial, prior = prior)
#' samples <- draw(sampler, n = 100, burnin = 10)
#'
#' x <- bnpostmcmc(sampler, samples)
#' ep(x)
priorGraph <- function(graph, lambda){
  stopifnot("bn" %in% class(graph),
            inherits(lambda, "numeric") || inherits(lambda, "integer"))
  lapply(seq_along(graph), function(i){
    i <- i # make the scope of i local
    function(x){
      difference <- setdiff(x, graph[[i]])
      difference <- length(difference)
      exp(- lambda * difference)
    }
  })
}

#' A uniform prior for graphs.
#'
#' A 'flat' improper prior that assigns equal probability to all the graphs.
#'
#' @param graph The 'prior graph'. A \code{bn}.
#' @return A function computes the prior score of the supplied graph. This
#'   This function is of a suitable form to be used as a prior
#' @export
#' @seealso \code{\link{priorGraph}}
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
priorUniform <- function(graph){
  stopifnot(inherits(graph, "bn"))
  lapply(seq_along(graph), function(x){
    function(x){
      1
    }
  })
}

#' Check validity.
#'
#' Checks the output of a LOG prior for validity. Basically is it positive?
#'
#' @param x A numeric of length 1. The LOG output of a prior function.
#' @return A logical of length 1, indicating the validity of x
#' @S3method is.valid prior
#' @method is.valid prior
#' @export
is.valid.prior <- function(x){
  tryCatch({
    all(class(x) %in% c("numeric", "integer"),
        !is.nan(x),
        is.finite(x))
    },
    error = function(e) F)
}


#' Check validity.
#'
#' Checks whether the supplied R object is a valid \code{localPriors}
#'
#' @param x An R object
#' @return A logical of length 1, indicating the validity of x
#' @export
is.valid.localPriors <- function(x){
  tryCatch({
    all(inherits(x, "list"),
        sapply(x, is.function))
    },
    error = function(e) F)
}

#' Compute prior score of a Bayesian network.
#'
#' Computes the complete prior of a network, from a \code{localPriors} object
#'
#' @param x An object of class 'bn'.
#' @param localPriors A list of functions of the same length as \code{x}
#'   that returns the local prior score of the corresponding node, given a
#'   numeric vector of parents.
#' @return The prior of the complete network \code{x}
#' @export
eval.prior <- function(x, localPriors){
  s <- seq_along(x)
  locals <- sapply(s, function(i){
    localPriors[[i]](x[[i]])
  })
  sum(locals)
}
