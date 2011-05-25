# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' BN Posterior from MCMC.
#'
#' method description
#'
#' @param sampler ...
#' @param samples ...
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
#' @seealso \code{\link{bnpost}}, \code{\link{bnpostmcmc.list}}.
#'   \code{\link{length.bnpostmcmc}}, \code{\link{top.bnpostmcmc}}, 
#'   \code{\link{map.bnpostmcmc}}, \code{\link{logScoreMultDir.bnpostmcmc}},
#'   \code{\link{gp.bnpostmcmc}}, \code{\link{ep.bnpostmcmc}}
bnpostmcmc <- function(sampler, samples, logScoreFUN){
  stopifnot(inherits(sampler, "sampler"),
            "bn.list"      %in% class(samples))

  if (get("return", envir = environment(sampler)) == "contingency"){
    tabulated <- get("count", envir = environment(sampler))
    tabulated <- as.list(tabulated)
    tabulated <- unlist(tabulated)
    tabulated <- sort.int(tabulated, method = "shell")
    tabulated <- as.table(tabulated)
  } else {
    tabulated <- pltabulate(samples, sort = T)
  }

  out <- new.env()
  out[["type"]] <- "mcmc"
  out[["logScoreFUN"]] <- evalq(expr  = logScoreFUN$offline,
                                envir = environment(sampler))
  out[["sampler"]] <- sampler
  out[["samples"]] <- samples
  out[["tabulated"]] <- tabulated
  out[["data"]] <- evalq(expr  = data,
                         envir = environment(sampler))

  class(out) <- "bnpostmcmc"
  out
}

#' List of BN Posteriors from MCMC.
#'
#' method description
#'
#' @param ... Further arguments passed to method
#' @S3method bnpostmcmc list
#' @method bnpostmcmc list
#' @export
#' @seealso \code{\link{bnpostmcmc}}, \code{\link{gp.bnpostmcmc.list}},
#'   \code{\link{ep.bnpostmcmc.list}}
bnpostmcmc.list <- function(...){
  x <- list(...)
  class(x) <- "bnpostmcmc.list"
  x
}

#' Number of samples drawn.
#' 
#' Returns the number of samples draw in the supplied object of class
#' 'bnpostmcmc' x.
#'
#' @param x An object of class 'bnpostmcmc'
#' @param ... Further arguments passed to method
#' @return The number of samples in x
#' @S3method length bnpostmcmc
#' @method length bnpostmcmc
#' @seealso \code{\link{bnpostmcmc}}.
length.bnpostmcmc <- function(x, ...){
  stopifnot(class(x) == "bnpostmcmc")
  length(x$samples)
}

#' Top graph from BN Posterior.
#' 
#' Returns the most commonly encountered graphs during the MCMC sampling
#' 'x'. The top 'head' graphs with respect to MCMC sampling are returned.
#' ie if the MCMC sampler has converged, the top graphs with respect to the
#' posterior distribution on graphs will be returned.
#'
#' @param x An object of class 'bnpostmcmc'
#' @param head The top head graphs will be returned. If the 'head'th
#'         most-commonly encountered graph ties in frequency with other
#'         graphs, all of the ties will be returned.
#' @param ... Further arguments (unused)
#' @return if head == 1:
#'     EITHER an object of class 'bn' containing the most commonly
#'     encountered graph,
#'     OR an object of class 'bn.list' containing a list of the equally most-
#'     commonly encountered graphs (if two or more were equally most-
#'     commonly) encountered
#'   if head > 1:
#'     an object of class 'bn.list' containing a list of the equally most-
#'     commonly encountered graphs
#' @S3method top bnpostmcmc
#' @method top bnpostmcmc
#' @seealso \code{\link{top}}
top.bnpostmcmc <- function(x, head = 10, ...){
  stopifnot(class(x) == "bnpostmcmc",
            is.wholenumber(head) || is.infinite(head),
            head > 0)
  tabulated <- x$tabulated
  rk <- rank(-tabulated, ties.method = "min") <= head
  top <- sort(tabulated[rk], dec = T)
  ids <- names(top)

  out <- as.parental(ids)
  if (length(ids) > 1){
    out <- lapply(out, function(bn){
      class(bn) <- c("bn", "parental")
      bn
    })
    class(out) <- c("bn.list", "parental.list")
  }
  else {
    class(out) <- c("bn", "parental")
  }
  out
}

#' Maximum aposteriori graph.
#' 
#' Returns the most commonly encountered graph(s) during the MCMC sampling
#' 'x'. ie the maximum aposteriori graph(s).
#'
#' @param x An object of class 'bnpostmcmc'
#' @param ... Further arguments (unused)
#' @return
#'   EITHER an object of class 'bn' containing the most commonly
#'   encountered graph,
#'   OR an object of class 'bn.list' containing a list of the equally most-
#'   commonly encountered graphs (if two or more were equally most-
#'   commonly) encountered
#' @S3method map bnpostmcmc
#' @method map bnpostmcmc
#' @seealso \code{\link{map}}
map.bnpostmcmc <- function(x, ...){
  stopifnot(class(x) == "bnpostmcmc")
  top(x, head = 1)
}

#' Log scores of best graphs.
#' 
#' Returns the log scores of the best graphs encountered by the sampler.
#'
#' @param x An object of class 'bnpostmcmc'
#' @param sampler A function generated by BNSampler() or similar
#' @param data A data.frame of the form required by logScoreMultDir()
#' @param sort.by A character of length 1. Either "posterior" or
#'              "logScoreMultDir".
#'              Determines the sorting of the result, and if head is finite
#'              and smaller than the number of graphs encountered in the
#'              MCMC, it also determines the graphs for which the scores are
#'              returned.
#'              if sort.by == "posterior":
#'                the scores of the most frequently encountered graphs
#'                are returned, sorted by the frequency that they are
#'                encountered. ie if the MCMC has converged, this returns
#'                best graphs with respect to the posterior distribution on
#'                graphs, sorted with respect to the posterior distribution
#'                on graphs.
#'              if sort.by == "logScoreMultDir":
#'                the scores of the graphs encountered during MCMC that
#'                score the most highly according to logScoreMultDir() are
#'                returned. This will be slower, because the scores of ALL
#'                the graphs encountered must be computed first (only local,
#'                incremental scores are stored by the MCMC)
#' @param head A numeric of length 1. The 1st to headth most-highly scoring
#'              (with respect to the sort.by metric) will be returned.
#'              If the 'head'th most-commonly encountered graph ties in
#'              score with other graphs, all of the ties will be returned.
#'              Thus the result may NOT be of length head.
#' @param use.names A logical of length 1. If TRUE, the result is a named vector,
#'              with names given by as.character(bn)
#' @param ... Further arguments passed to logScoreMultDir()
#'
#' @return A vector of length at least head (if head is finite), giving the
#'   logScoreMultDir() of the graphs, or all the graphs encountered during
#'   MCMC if head == Inf.
#' @S3method logScoreMultDir bnpostmcmc
#' @method logScoreMultDir bnpostmcmc
#' @seealso \code{\link{logScoreMultDir}}
logScoreMultDir.bnpostmcmc <- function(x, sampler, data,
                                       sort.by   = "posterior",
                                       head      = Inf,
                                       use.names = F, ...){
  stopifnot(class(x)          ==   "bnpostmcmc",
            inherits(sampler, "sampler"),
            class(data)       ==   "data.frame",
            is.wholenumber(head) || is.infinite(head),
            head              >    0,
            sort.by           %in% c("posterior", "logScoreMultDir"),
            class(use.names)  ==   "logical",
            length(use.names) ==   1)
  # extract the cache from the sampler to speed up computation
  cache <- get("cache", envir = environment(sampler))
  # we must score all the graphs if we sort by logScoreMultDir
  if (sort.by == "logScoreMultDir"){
    bnptop <- top(x, head = Inf)
  }
  else {
    bnptop <- top(x, head = head)
  }
  res <- logScoreMultDir(bnptop, data, cache = cache, ...)
  if ("bn" %in% class(bnptop)){
    # TODO: this seems like a bit of a fiddle.
    res <- res[[1]]
  }
  if (isTRUE(use.names)){
    names(res) <- as.character(bnptop)
  }
  if (sort.by == "logScoreMultDir"){
    rk <- rank(-res, ties.method = "min") <= head
    res <- sort(res[rk], dec = T)
  }
  res
}

#' Posterior graph probabilities.
#'
#' method description
#'
#' @param x ...
#' @param start ...
#' @param end ...
#' @param nbin ...
#' @param log ...
#' @param pretty ...
#' @param levels ...
#' @param ... Further arguments (unused)
#' @S3method gp bnpostmcmc
#' @method gp bnpostmcmc
#' @seealso \code{\link{gp}}, \code{\link{gp.bnpostmcmc.list}}
gp.bnpostmcmc <- function(x, start, end, nbin = 1,
                          log = F, pretty = F, levels = NULL, ...){

  makeNames <- function(tabulatedProportions){
    # a probably suboptimal shuffle so as to keep names
    if (pretty){
      bnl <- as.parental(names(tabulatedProportions))
      ids <- as.character(bnl, pretty = T)
    }
    else {
      ids <- names(tabulatedProportions)
    }
    tabulatedProportions <- as.vector(tabulatedProportions)
    names(tabulatedProportions) <- ids
    tabulatedProportions
  }

  if (nbin == 1){
    nSamples <- length(x$samples)
    out <- x$tabulated/nSamples

    out <- makeNames(out)

    class(out) <- "gp"
    out
  }
  else {
    # if the nbin > 1, we have to retabulate
    lengthOfRuns <- length(x$samples)
    sizeOfBins <- lengthOfRuns/nbin
    if (round(sizeOfBins) != sizeOfBins) stop("nbin not an integer malarkey")
    samples <- x$samples

    #bins <- split(samples, rep(seq_len(nbin), each = sizeOfBins))
    if (is.null(levels)){
      levels <- names(x$tabulated)
    }

    # loop over the bins
    outs <- lapply(seq_len(nbin), function(bin){
      samples <- samples[
        seq.int(
          from = sizeOfBins * (bin - 1) + 1,
          to = sizeOfBins * bin
        )
      ]
      class(samples) <- c("parental.list")

      # need to keep this sorted according to the levels
      out <- pltabulate(samples, levels = levels,
                        pretty = pretty, sort = F)
      out <- out/length(samples)

      #out <- makeNames(out, pretty)
      class(out) <- "gp"
      out
    })
    class(outs) <- "gp.list"
    outs
  }
}

#' Posterior graph probabilities.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @S3method gp bnpostmcmc.list
#' @method gp bnpostmcmc.list
#' @seealso \code{\link{gp}}, \code{\link{gp.bnpostmcmc}}
gp.bnpostmcmc.list <- function(x, ...){
  stopifnot(class(x) == "bnpostmcmc.list")
  lapply(x, gp, ...)
}

#' Posterior edge probabiities.
#' 
#' Computes the edge probabilities implied by the MCMC samples
#' contained in the 'bnpostmcmc' object x.
#'
#' @param x An object of class 'bnpostmcmc'
#' @param nbin A numberic vector of length 1 specifying the number of bins
#'           into which to divide the MCMC samples. The edge probabilities
#'           are computed separately for each bin.
#' @param start ...
#' @param end ...
#' @param method Either "et" (the default), "flatten" or "tabulate". If
#'   the edge totals are not available, method "tabulate" is used if
#'   possible. Only "flatten" is available if \code{nbin != 1}.
#' @param verbose ...
#' @param ... Further arguments passed to ep.parental.list() for method =
#'           "flatten", or ep.table() for method = "tabulate"
#' @return if nbin == 1:
#'     A matrix of class 'ep' with entry (i,j) containing the probability of
#'     an edge from node i --> j
#'   if nbin > 1:
#'     A list of class ep.list, containing matrices as described above for
#'     each of the nbin bins into which the parental.list was split
#' @S3method ep bnpostmcmc
#' @method ep bnpostmcmc
#' @seealso \code{\link{ep}}, \code{\link{ep.bnpostmcmc.list}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#' 
#' prior <- function(net) 1
#' initial <- bn(c(), c())
#' 
#' sampler <- BNSampler(dat, initial, prior)
#' samples <- draw(sampler, n = 50)
#' mpost <- bnpostmcmc(sampler, samples)
#' 
#' ep(mpost)
ep.bnpostmcmc <- function(x, nbin = 1, start, end, method = "et",
                          verbose = F, ...){
  stopifnot(class(x) == "bnpostmcmc",
            isTRUE(is.wholenumber(nbin)),
            method %in% c("et", "flatten", "tabulate"))
  nbin_1 <- nbin == 1
  etbins_exists <- exists("etbins", envir = environment(x$sampler))
  if (method != "flatten" && !nbin_1){
    stop("Only method = 'flatten' is implemented for nbin != 1")
  }

  if (method == "et" && nbin_1 && etbins_exists){
    if (verbose){
      cat("Using edge total matrix to compute ep\n")
    }
    et <- get("etbins", envir = environment(x$sampler))
    numberOfNodes <- get("numberOfNodes", envir = environment(x$sampler))
    et <- matrix(colSums(et, na.rm = T), numberOfNodes, numberOfNodes)
    nSteps <- get("nSteps", envir = environment(x$sampler))
    nBurnin <- get("nBurnin", envir = environment(x$sampler))
    ep <- et/(nSteps - nBurnin)
    class(ep) <- c("ep", "matrix")
    ep
  } else if (method == "tabulate" && nbin_1){
    if (verbose){
      cat("Using tabulated samples to compute ep\n")
    }
    ep(x$tabulated, verbose = verbose, ...)
  }  else if (method == "flatten"){
    if (verbose){
      cat("Flattening the samples to compute ep\n")
    }
    ep(x$samples, nbin, ...)
  }
}

#' Extract posterior edge probabiities from a sampler.
#' 
#' Computes the edge probabilities implied by the MCMC samples
#' contained in a MCMC sampler function.
#'
#' @param x An MCMC sampler 'function'
#' @param start ...
#' @param end ...
#' @param verbose ...
#' @param ... Further arguments passed to ep.parental.list() for method =
#'           "flatten", or ep.table() for method = "tabulate"
#' @return A matrix of class 'ep' with entry (i,j) containing the probability
#'   of an edge from node i --> j
#' @S3method ep sampler
#' @method ep sampler
#' @seealso \code{\link{ep}}, \code{\link{ep.bnpostmcmc.list}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#'
#' prior <- function(net) 1
#' initial <- bn(c(), c())
#'
#' sampler <- BNSampler(dat, initial, prior)
#' samples <- draw(sampler, n = 50)
#'
#' ep(sampler)
ep.sampler <- function(x, start, end, verbose = F, ...){
  stopifnot(inherits(x, "sampler"))

  etbins_exists <- exists("etbins", envir = environment(x))

  if (etbins_exists){
    if (verbose){
      cat("Using edge total matrix to compute ep\n")
    }
    et <- get("etbins", envir = environment(x))
    numberOfNodes <- get("numberOfNodes", envir = environment(x))
    et <- matrix(colSums(et, na.rm = T), numberOfNodes, numberOfNodes)
    nSteps <- get("nSteps", envir = environment(x))
    nBurnin <- get("nBurnin", envir = environment(x$sampler))
    ep <- et/(nSteps - nBurnin)
    class(ep) <- c("ep", "matrix")
    ep
  } else {
    stop("The supplied function 'x' does not appear to be a sampler.")
  }
}

#' Posterior edge probabilities.
#'
#' method description
#'
#' @param x ...
#' @param start ...
#' @param end ...
#' @param ... further arguments
#' @S3method ep bnpostmcmc.list
#' @method ep bnpostmcmc.list
#' @seealso \code{\link{ep}}, \code{\link{ep.bnpostmcmc}}
#' @examples
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#' 
#' prior <- function(net) 1
#' initial <- bn(c(), c())
#' sampler <- BNSampler(dat, initial, prior)
#' samples <- draw(sampler, n = 50)
#' mpost <- bnpostmcmc(sampler, samples)
#' 
#' initial <- bn(c(), c(1))
#' sampler2 <- BNSampler(dat, initial, prior)
#' samples2 <- draw(sampler2, n = 50)
#' mpost2 <- bnpostmcmc(sampler2, samples2)
#' 
#' ep(bnpostmcmc.list(mpost, mpost2))
ep.bnpostmcmc.list <- function(x, start, end, ...){
  stopifnot(class(x) == "bnpostmcmc.list")
  if (!missing(start) || !missing(end)){
    if (start != 1 && end != length(x[[1]])){
      warning("Start/end not implemented")
    }
  }
  # FIXME: using ep.bnpostmcmc is a temporary fix
  # otherwise method dispatches to ep.parental.list
  res <- lapply(x, function(y){
    ep(y, ...)
  })
  class(res) <- "ep.list"
  res
}
