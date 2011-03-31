# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Posterior distribution on Bayesian networks
#'
#' Use one of a number of methods to get the posterior distribution.
#'
#' @param data The data.
#' @param method One of "exact", "mh-mcmc", "gibbs", "mj-mcmc"
#' @param prior A function that returns the prior score of the
#'                       supplied bn.
#' @param logScoreFUN A list of three elements:
#'                         offline: A function that computes the logScore
#'                                  of a Bayesian Network
#'                         online:  A function that incrementally computes
#'                                  the logScore of a Bayesian Network
#'                         prepare: A function that prepares the data, and
#'                                  any further pre-computation required by
#'                                  the logScore functions.
#' @param logScoreParameters A list of parameters that are passed to
#'                       logScoreFUN.
#' @param constraint A matrix of dimension ncol(data) x ncol(data) giving
#'                       constraints to the sample space.
#'                       The (i, j) element is
#'                         1  if the edge i -> j is required
#'                         -1 if the edge i -> is excluded.
#'                         0  if the edge i -> j is not constrained.
#'                       The diagonal of constraint must be all 0.
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node.
#' @param nSamples The number of samples to be draw (only applies to MCMC)
#' @param nBurnin The number of samples to discard from the beginning of
#'   the sample.
#' @param initial An object of class 'bn'. The starting value of the
#'                       MCMC.
#' @param verbose A logical. Should a progress bar be displayed?
#' @return Either a \code{bnpost} or a \code{bnpostmcmc} object.
#' @export
posterior <- function(data,
                      method             = "mh-mcmc",
                      prior              = function(x) 1,
                      logScoreFUN        = defaultLogScoreFUN(),
                      logScoreParameters = list(hyperparameters = "qi"),
                      constraint         = NULL,
                      maxNumberParents   = Inf,
                      nSamples           = 50000,
                      nBurnin            = 10000,
                      initial            = NULL,
                      verbose            = F){
  methods <- c("exact", "mh-mcmc", "gibbs", "mj-mcmc")
  stopifnot(class(data) ==   "data.frame",
            method      %in% methods)

  if (method == "exact"){
    exactposterior(data,
                   prior,
                   logScoreFUN,
                   logScoreParameters,
                   constraint,
                   maxNumberParents,
                   verbose)
  } else if (method == "mh-mcmc"){
    mcmcposterior(sampler = BNSampler,
                  data,
                  prior,
                  logScoreFUN,
                  logScoreParameters,
                  constraint,
                  maxNumberParents,
                  nSamples,
                  nBurnin,
                  initial,
                  verbose)
  } else if (method == "gibbs") {
    mcmcposterior(sampler = BNGibbsSampler,
                  data,
                  prior,
                  logScoreFUN,
                  logScoreParameters,
                  constraint,
                  maxNumberParents,
                  nSamples,
                  nBurnin,
                  initial,
                  verbose)
  } else if (method == "mj-mcmc") {
    mcmcposterior(sampler = BNSamplerMJ,
                  data,
                  prior,
                  logScoreFUN,
                  logScoreParameters,
                  constraint,
                  maxNumberParents,
                  nSamples,
                  nBurnin,
                  initial,
                  verbose)
  } else {
    stop("Not implemented")
  }
}

#' Posterior distribution on Bayesian networks
#'
#' Use one of a number of methods to get the posterior distribution
#'
#' @param data The data.
#' @param prior A function that returns the prior score of the
#'                       supplied bn.
#' @param logScoreFUN A list of three elements:
#'                         offline: A function that computes the logScore
#'                                  of a Bayesian Network
#'                         online:  A function that incrementally computes
#'                                  the logScore of a Bayesian Network
#'                         prepare: A function that prepares the data, and
#'                                  any further pre-computation required by
#'                                  the logScore functions.
#' @param logScoreParameters A list of parameters that are passed to
#'                       logScoreFUN.
#' @param constraint A matrix of dimension ncol(data) x ncol(data) giving
#'                       constraints to the sample space.
#'                       The (i, j) element is
#'                         1  if the edge i -> j is required
#'                         -1 if the edge i -> is excluded.
#'                         0  if the edge i -> j is not constrained.
#'                       The diagonal of constraint must be all 0.
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node.
#' @param verbose A logical. Should a progress bar be displayed?
#' @return A \code{bnpost} object.
#' @export
exactposterior <- function(data,
                           prior,
                           logScoreFUN,
                           logScoreParameters,
                           constraint,
                           maxNumberParents,
                           verbose){
  nVar <- ncol(data)
  bnspace <- enumerateBNSpace(nVar)

  if (!is.null(constraint)){
    bnspace <- Filter(satisfiesConstraint, bnspace)
  }

  logScoreOfflineFUN <- logScoreFUN$offline
  prepareDataFUN <- logScoreFUN$prepare
  logScoreParameters <- prepareDataFUN(data,
                                       logScoreParameters,
                                       checkInput = F)

  if (isTRUE(verbose)){
    progress <- create_progress_bar("text")
    progress$init(length(bnspace))
  }
  
  lsmd <- sapply(bnspace, function(x){
    if (isTRUE(verbose)){
      progress$step()
    }
    logScoreOfflineFUN(x                  = x,
                       logScoreParameters = logScoreParameters)
  })
  
  if (isTRUE(verbose)){
    progress$term()
  }

  logpriors <- log(sapply(bnspace, prior))
  logScore <- lsmd + logpriors
  
  bnpost(bnspace     = bnspace,
         logScore    = logScore,
         data        = data,
         logScoreFUN = function(x) stop("error"))
}

#' Posterior distribution on Bayesian networks
#'
#' Use MCMC to approximate the posterior distribution
#'
#' @param sampler A BNSampler. eg BNSampler or BNGibbsSampler etc
#' @param data The data.
#' @param prior A function that returns the prior score of the
#'                       supplied bn.
#' @param logScoreFUN A list of three elements:
#'                         offline: A function that computes the logScore
#'                                  of a Bayesian Network
#'                         online:  A function that incrementally computes
#'                                  the logScore of a Bayesian Network
#'                         prepare: A function that prepares the data, and
#'                                  any further pre-computation required by
#'                                  the logScore functions.
#' @param logScoreParameters A list of parameters that are passed to
#'                       logScoreFUN.
#' @param constraint A matrix of dimension ncol(data) x ncol(data) giving
#'                       constraints to the sample space.
#'                       The (i, j) element is
#'                         1  if the edge i -> j is required
#'                         -1 if the edge i -> is excluded.
#'                         0  if the edge i -> j is not constrained.
#'                       The diagonal of constraint must be all 0.
#' @param maxNumberParents Integer of length 1. The maximum number of
#'   parents of any node.
#' @param nSamples The number of samples to be draw (only applies to MCMC)
#' @param nBurnin The number of samples to discard from the beginning of
#'   the sample.
#' @param initial An object of class 'bn'. The starting value of the
#'                       MCMC.
#' @param verbose A logical. Should a progress bar be displayed?
#' @return A \code{bnpostmcmc} object.
#' @export
mcmcposterior <- function(sampler = BNSampler,
                          data,
                          prior,
                          logScoreFUN,
                          logScoreParameters,
                          constraint,
                          maxNumberParents,
                          nSamples,
                          nBurnin,
                          initial,
                          verbose){
  nVar <- ncol(data)
  if (is.null(initial)){
    initial <- empty(nVar, class = "bn")
  }

  sampler <- sampler(data               = data,
                     initial            = initial,
                     prior              = prior,
                     logScoreFUN        = logScoreFUN,
                     logScoreParameters = logScoreParameters,
                     constraint         = constraint,
                     maxNumberParents   = maxNumberParents,
                     verbose            = verbose)
  samples <- draw(sampler = sampler,
                  n       = nSamples,
                  burnin  = nBurnin,
                  verbose = verbose)
  bnpostmcmc(sampler     = sampler,
             samples     = samples,
             logScoreFUN = logScoreFUN)
}

#' method name
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
gp <- function(x, ...){
  UseMethod("gp")
}

#' method name
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
ep <- function(x, ...){
  UseMethod("ep")
}

#' method name
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
entropy <- function(x, ...){
  UseMethod("entropy")
}

#' method name
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
map <- function(x, ...){
  UseMethod("map")
}

#' method name
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
top <- function(x, ...){
  UseMethod("top")
}

#' method name
#'
#' method description
#'
#' @param bn1 A bn
#' @param bn2 A bn
#' @param data data
#' @param ... Further arguments passed to method
#' @export
bf <- function(bn1, bn2, data, ...){
  UseMethod("bf")
}

#' Calculate edge probabilities given a list of
#' parental graphs (of class parental.list)
#'
#' @param x A parental.list
#' @param nbin The number of equally-sized bins into which parental.list into
#'          before computing the edge probabilities of each
#' @param start An integer of length 1, specifying the amount of burn-in.
#'          The samples start:end inclusive will be used.
#' @param end An integer of length 1, specifying the number of samples at the
#'          end to ignore. The samples start:end inclusive will be used.
#' @param verbose ...
#' @param ... Further arguments (unused)
#' @return if nbin == 1:
#'     A matrix of class 'ep' with entry (i,j) containing the probability of
#'     an edge from node i --> j
#'   if nbin > 1:
#'     A list of class ep.list, containing matrices as described above for
#'     each of the nbin bins into which the parental.list was split
#' @S3method ep parental.list
#' @method ep parental.list
ep.parental.list <- function(x, nbin = 1, start = 1, end = length(x),
                             verbose = F, ...){
  stopifnot("parental.list" %in% class(x),
            isTRUE(is.wholenumber(nbin)),
            validStartEnd(start, end, length(x)))
  # remove burn-in etc
  x <- x[seq(from = start, to = end)]

  lengthOfRun <- length(x)
  binSeq <- seq_len(nbin)
  numberOfNodes <- length(x[[1]])
  nodeSeq <- seq_len(numberOfNodes)
  sizeOfBins <- lengthOfRun/nbin
  if (!isTRUE(is.wholenumber(sizeOfBins))){
    stop("The number of bins (nbin) must divide into the number of samples")
  }

  # flatten samples early for efficiency
  # and remove samples for memory efficiency
  if (verbose){
    cat("Flattening samples\n")
  }

  #flatsamples <- vector("list", length = lengthOfRun * numberOfNodes);
  x <- unlist(x, recursive = F, use.names = F);
  #x <- NULL

  extractSamples <- function(bin){
    # By flattening the data and picking out every numberOfNodes th
    # element, we get the appropriate vector that shows the parents of
    # head.

    # create a subsetting vector:
    # We have data of the form
    #
    # list(
    #   list(2, integer(0), 1),
    #   list(integer(0), 1, 1),
    #   list(integer(0), 1, 1),
    #   list(integer(0), integer(0), 1),
    #   list(integer(0), 1, 1),
    #   list(integer(0), integer(0), 1)
    # )
    #
    # if nbin = 3
    # then for bin = 2, head = 1
    # we want
    #
    # from = 2 * 3 * (2 - 1) + 1 = 6 + 1  = 7
    # to = 2 * 3 * (2 - 1) + 1 + 2 * 3 - 1 = 6 + 1 + 6 - 1 = 12
    # which gives
    # 7, 8, 9, 10, 11, 12
    # as required

    tabulateNULL <- function(bin, ...){
      if (is.null(bin)){
        rep.int(0, numberOfNodes)
      }
      else {
        tabulate(bin, ...)
      }
    }

    out <- matrix(ncol = numberOfNodes, nrow = numberOfNodes)
    toTabulate <- vector("numeric", length = lengthOfRun * numberOfNodes)

    if (verbose){
      cat("Iterating over nodes\n")
    }
    for (head in nodeSeq){
      base <- sizeOfBins * numberOfNodes * (bin - 1) + head
      select <- seq.int(
        from = base,
        to = base + (sizeOfBins - 1) * numberOfNodes,
        by = numberOfNodes
      )
      if (verbose){
        cat("Unlisting data\n")
      }
      toTabulate <- unlist(x[select],
                           recursive = F,
                           use.names = F)

      # use tabulateNULL because if head never has any parents in this bin, 
      # then flatten will contain a NULL and tabulate throws an error on NULL.
      # Place the output in every numberOfNodes th column
      if (verbose){
        cat("Tabulating\n")
      }
      out[, head] <- tabulateNULL(toTabulate, nbins = numberOfNodes)
    }
    epmx <- out/sizeOfBins
    class(epmx) <- c("ep", "matrix")
    epmx
  }
  if (verbose){
    cat("Extracting samples\n")
  }
  epl <- lapply(binSeq, extractSamples)

  if (nbin == 1){
    eparr <- epl[[1]]
    class(eparr) <- c("ep", "matrix")
    eparr
  }
  else {
    class(epl) <- "ep.list"
    epl
  }
}

#' Computes the edge probabilities implied by a table.
#'
#' Computes the edge probabilities implied by a table.
#' 
#' @param x An object of class 'table'. The table should have names that
#'   can be converted into objects of class parental using
#'   \code{as.parental.character}. ie they should be serializations of
#'   parental objects.
#' @param verbose A logical indicating whether progress should be shown.
#' @param ... Further arguments passed to
#'   \code{\link{ep.parental.contingency}}
#' @return A matrix of class 'ep' with entry (i,j) containing the probability
#'   of an edge from node i --> j
#' @S3method ep table
#' @method ep table
ep.table <- function(x, verbose = F, ...){
  stopifnot(class(x) == "table")

  tablenames <- names(x)
  if (verbose) cat("Converting names to parental\n")
  tablebn <- as.parental(tablenames)
  
  pc <- list(parental.list = tablebn,
             contingency   = as.numeric(x))
  class(pc) <- "parental.contingency"
  ep(pc, verbose = verbose, ...)
}

#' Computes the edge probabilities
#' 
#' Computes the edge probabilities from a \code{parental.contingency},
#' which is a list, whose first component is \code{parental.list}, and
#' whose second component is a contingency table of counts of the
#' corresponding element of the \code{parental.list}. The contingency
#' table must simply be a numeric vector. The first component must be named
#' \code{parental.list}, and the second \code{contingency}.
#'
#' @param x An object of class \code{parental.contingency}, of the form
#'   desribed in the Details section.
#' @param FUN A function to be applied to the samples before the edge
#'   probabilities are computed. The function should accept a vector
#'   of objects of class 'parental.list' and return the transformed
#'   parental objects as a single 'parental.list'. If the function
#'   does not return a 'parental.list', an error will be thrown. Can be
#'   supplied as a function, symbol or character string (see
#'   \code{\link[base]{match.fun}})
#' @param verbose A logical indicating whether progress should be shown.
#' @param ... Further arguments (unused)
#' @return A matrix of class 'ep' with entry (i,j) containing the
#'  probability of an edge from node \code{i} to \code{j}
#' @S3method ep parental.contingency
#' @method ep parental.contingency
ep.parental.contingency <- function(x, FUN, verbose = F, ...){
  stopifnot(class(x) == "parental.contingency",
            "parental.list" %in% names(x),
            "contingency" %in% names(x),
            inherits(x$parental.list, "parental.list"),
            inherits(x$contingency, c("numeric", "integer")),
            inherits(verbose, "logical"))

  tablebn <- x$parental.list
  counts <- x$contingency
  
  tablebn <- unname(tablebn)
  counts <- unname(counts)
  
  # apply FUN if supplied
  if (!missing(FUN)){
    if (verbose) cat("Applying FUN to parental.list\n")
    FUN <- match.fun(FUN)
    tablebn <- FUN(tablebn)
    stopifnot(class(tablebn) == "parental.list")
  }

  numberOfNodes <- length(tablebn[[1]])
  nodeSeq <- seq_len(numberOfNodes)
  lengthOfRuns <- sum(counts)

  if (verbose) cat("Computing ep() from parental.contingency\n")
  res <- matrix(0, ncol = numberOfNodes, nrow = numberOfNodes)
  # loop over each cell of the edge probabilities matrix
  if (verbose){
    progress <- create_progress_bar("text")
    progress$init(numberOfNodes^2)
  }
  
  # t <- system.time({
  # for (i in seq_along(tablebn)){
  #   for (head in nodeSeq){
  #     tails <- tablebn[[i]][[head]]
  #     res[tails, head] <- res[tails, head] + counts[i]
  #   }
  #   # if (verbose){
  #   #   progress$step()
  #   # }
  # }})
  # 
  # browser()

  for (head in nodeSeq){
    thishead <- lapply(tablebn, "[[", head)
    lenhead <- sapply(thishead, length)
    counthead <- rep(counts, times = lenhead)
    thishead <- unlist(thishead)
    for (tail in nodeSeq){
      res[tail, head] <- sum(counthead[thishead == tail])
      if (verbose){
        progress$step()
      }
    }
  }
  if (verbose) cat("\n")
  
  res <- res/lengthOfRuns
  class(res) <- c("ep", "matrix")
  res
}

#' Plot the the top() Bayesian Networks from a posterior distribution.
#' The top graphs are those with the highest score with respect to the
#' posterior distribution, which for converged MCMC will be most commonly
#' encountered graph.
#'
#' @param x An object of class "bnpostmcmc" or "bnpost"
#' @param top Optionally provide pre-computed top(x)
#'
#' @return A plot of the top graph, with their marginal likelihoods (without
#'   priors)
#' @S3method plot bnpostmcmc
#' @method plot bnpostmcmc
plot.bnpostmcmc <- plot.bnpost <- function(x, top = NULL, ...){
  if (is.null(top)){
    top <- top(x)
  }
  lsmd <- logScoreMultDir(top, data = x$data)
  lsmd <- round(lsmd, 2)
  names(top) <- paste(seq_along(top), ": ", as.character(lsmd), sep = "")
  grplot(top, ...)
}

#' method name
#'
#' method description
#'
#' @param bn1 A bn
#' @param bn2 ...
#' @param data ...
#' @param ... further arguments
#'
#' @S3method bf bnpostmcmc
#' @method bf bnpostmcmc
bf.bnpostmcmc <- bf.bnpost <- function(bn1, bn2, data, ...){
  bnl <- bn.list(bn1, bn2)
  logScoreMultDir(bnl, data, ...)
}

#' method name
#'
#' method description
#'
#' @param x ...
#'
#' @S3method levelplot bnpostmcmc
#' @method levelplot bnpostmcmc
levelplot.bnpostmcmc <- levelplot.bnpost <- function(x){
  stopifnot(class(x) %in% c("bnpostmcmc", "bnpost"))
  levelplot(ep(x))
}

#' method name
#'
#' method description
#'
#' @param ep ...
prepareLevelPlot <- function(ep){
  n <- dim(ep)[1]
  nodeSeq <- seq_len(n)

  data <- expand.grid(head = factor(nodeSeq, levels = rev(nodeSeq)),
                      tail = factor(nodeSeq, levels = nodeSeq))

  data$edgeprob <- as.vector(as.numeric(ep))

  data
}

#' method name
#'
#' method description
#'
#' @param ep ...
#'
#' @S3method levelplot ep
#' @method levelplot ep
levelplot.ep <- function(ep){
  data <- prepareLevelPlot(ep)

  levelplot(
    edgeprob ~ tail * head,
    data,
    ylab = "Tail",
    xlab = "Head",
    main = "Edge probabilities",
    panel = function(x, y, z, subscripts, ...){
      panel.levelplot(x, y, z, subscripts, ...)
      ltext(x, y, signif(z[subscripts], 2), col = "black")
    }
  )
}

#' method name
#'
#' method description
#'
#' @param eplist ...
#' @S3method levelplot ep.list
#' @method levelplot ep.list
levelplot.ep.list <- function(eplist){
  dfs <- lapply(eplist, prepareLevelPlot)
  data <- do.call("make.groups", dfs)

  levelplot(
    edgeprob ~ tail * head | which,
    data,
    ylab = "Tail",
    xlab = "Head",
    aspect = "iso",
    as.table = T,
    main = "Edge probabilities",
    panel = function(x, y, z, subscripts, ...){
      panel.levelplot(x, y, z, subscripts, ...)
      ltext(x[subscripts], y[subscripts],
            signif(z[subscripts], 2), col = "black")
    }
  )
}

#' method name
#'
#' method description
#'
#' @param ep ...
shrinkep <- function(ep){
  offDiagonals <- upper.tri(ep) + lower.tri(ep) > 0

  # remove the diagonals
  ep <- ep[offDiagonals]

  # name the edge probabilities
  whichEdges <- which(offDiagonals, arr.ind = T)
  edgeNames <- apply(whichEdges, 1, paste, collapse = "->")

  names(ep) <- edgeNames
  ep
}

#' method name
#'
#' method description
#'
#' @param ep ...
#' @param head ...
#' @param ... Further arguments passed to method
#' @S3method dotplot ep
#' @method dotplot ep
dotplot.ep <- function(ep, head = 30, ...){
  epdata <- data.frame(method = c(), ep = c(), name = c())

  epShrunk <- shrinkep(ep)

  epdata <- rbind(epdata, data.frame(
    ep = epShrunk,
    name = names(epShrunk)
  ))

  subsetTotals <- tapply(epdata[, c("ep")], epdata[, "name"], sum)
  subsetNames <- names(sort(subsetTotals, dec = T))

  # this was reorder(name) before R-check picked up on problem
  epdata$name <- with(epdata, reorder(epdata$name, ep))

  if (!is.null(head)){
    subsetNames <- subsetNames[seq_len(head)]
  }

  dotplot(
    name ~ ep,
    data = epdata,
    type = c("p", "g"),
    xlab = "Edge probability",
    auto.key = T,
    ylab = "Edges",
    subset = epdata$name %in% subsetNames,
    par.settings = simpleTheme(pch = 3:10),
    ...
  )
}

#' method name
#'
#' method description
#'
#' @param ... Further arguments passed to method
#' @export
ep.list <- function(...){
  out <- list(...)
  class(out) <- "ep.list"
  out
}

#' method name
#'
#' method description
#'
#' @param eplist ...
#' @param subsetBy ...
#' @param head ...
#' @param ... Further arguments passed to method
#' @S3method dotplot ep.list
#' @method dotplot ep.list
dotplot.ep.list <- function(eplist, subsetBy = "Exact", head = 30, ...){
  epTypes <- names(eplist)

  epdataList <- lapply(seq_along(eplist), function(whichEP){
    ep <- eplist[[whichEP]]
    nm <- epTypes[[whichEP]]

    epShrunk <- shrinkep(ep)
    data.frame(
      method = nm,
      ep = epShrunk,
      name = names(epShrunk)
    )
  })

  epdata <- do.call("rbind", epdataList)

  epdata[, "method"] <- as.factor(epdata[, "method"])

  if (!subsetBy %in% epTypes){
    subsetBy <- "MCMC"
  }

  epdatasubset <- subset(epdata, epdata$method %in% subsetBy)
  subsetTotals <- tapply(epdatasubset[, c("ep")], epdatasubset[, "name"], sum)
  subsetNames <- names(sort(subsetTotals, dec = T))

  epdata$name <- with(epdata, reorder(epdata$name, ep))

  if (!is.null(head)){
    subsetNames <- subsetNames[seq_len(head)]
  }

  dotplot(
    name ~ ep,
    groups       = epdata$method,
    data         = epdata,
    type         = c("p", "g"),
    xlab         = "Edge probability",
    auto.key     = T,
    ylab         = "Edges",
    subset       = epdata$name %in% subsetNames,
    par.settings = simpleTheme(pch = 3:10),
    ...
  )
}

#' method name
#'
#' method description
#'
#' @param gplist ...
prepareGPPlot <- function(gplist){
  if (class(gplist) == "gp.list"){

    gpdataList <- lapply(names(gplist), function(nm){
      gp <- gplist[[nm]]

      data.frame(
        method = nm,
        gp = as.vector(gp),
        name = names(gp)
      )
    })
    gpdata <- do.call("rbind", gpdataList)

    gpdata[, "method"] <- as.factor(gpdata[, "method"])
    gpdata
  }
  else {
    # otherwise gplist is not a list

    data.frame(
        gp = as.vector(gplist),
        name = names(gplist)
    )
  }
}

#' method name
#'
#' method description
#'
#' @param gp ...
#' @param head ...
#' @param ... Further arguments passed to method
#' @S3method dotplot gp
#' @method dotplot gp
dotplot.gp <- function(gp, head = 30, ...){
  gpdata <- prepareGPPlot(gp)

  subsetTotals <- tapply(gpdata[, c("gp")], gpdata[, "name"], sum)
  subsetNames <- names(sort(subsetTotals, dec = T))

  gpdata$name <- with(gpdata, reorder(epdata$name, gp))

  if (!is.null(head)){
    subsetNames <- subsetNames[seq_len(head)]
  }

  dotplot(
    name ~ gp,
    data         = gpdata,
    type         = c("p", "g"),
    xlab         = "Graph probability",
    auto.key     = T,
    ylab         = "Graphs",
    subset       = gpdata$name %in% subsetNames,
    par.settings = simpleTheme(pch = 3:10),
    ...
  )
}

#' method name
#'
#' method description
#'
#' @param gp ...
#' @param head ...
#' @param ... Further arguments passed to method
#' @S3method xyplot gp
#' @method xyplot gp
xyplot.gp <- function(gp, head = 30, ...){
  gpdata <- prepareGPPlot(gp)

  # compute the probability of each graph
  subsetTotals <- tapply(gpdata[, c("gp")], gpdata[, "name"], sum)
  # sort the names by decreasing order of probability
  subsetNames <- names(sort(subsetTotals, dec = T))
  # sort the name factor by graph probability
  gpdata$name <- with(gpdata, reorder(name, gp))

  # reverse order of the levels
  gpdata$name <- factor(as.character(gpdata$name),
                        levels = rev(levels(gpdata$name)))

  if (!is.null(head)){
    subsetNames <- subsetNames[seq_len(head)]
  }

  xyplot(
    gp ~ name,
    data         = gpdata,
    type         = "h",
    xlab         = "Graphs",
    ylab         = "Graph probability",
    scales       = list(x = list(rot = 90)),
    subset       = gpdata$name %in% subsetNames,
    par.settings = simpleTheme(pch = 3:10),
    ...
  )
}

#' method name
#'
#' method description
#'
#' @param ... Further arguments passed to method
#' @export
gp.list <- function(...){
  out <- list(...)
  class(out) <- "gp.list"
  out
}

#' method name
#'
#' method description
#'
#' @param gplist ...
#' @param subsetBy ...
#' @param head ...
#' @param ... Further arguments passed to method
#' @S3method dotplot gp.list
#' @method dotplot gp.list
dotplot.gp.list <- function(gplist, subsetBy = "Exact", head = 30, ...){
  gpdata <- prepareGPPlot(gplist)

  if (!"Exact" %in% gpdata[, "method"]){
    subsetBy <- names(gplist)[1]
    warning(paste("Subsetting by", subsetBy))
  }

  gpdatasubset <- subset(gpdata, gpdata$method %in% subsetBy)
  subsetTotals <- tapply(gpdatasubset[, c("gp")], gpdatasubset[, "name"], sum)
  subsetNames <- names(sort(subsetTotals, dec = T))

  gpdata$name <- with(gpdata, reorder(name, gp))

  if (!is.null(head)){
    subsetNames <- subsetNames[seq_len(head)]
  }

  dotplot(
    name ~ gp,
    groups       = gpdata$method,
    data         = gpdata,
    type         = c("p", "g"),
    xlab         = "Graph probability",
    auto.key     = T,
    ylab         = "Graphs",
    subset       = gpdata$name %in% subsetNames,
    par.settings = simpleTheme(pch = 3:10),
    ...
  )
}

#' method name
#'
#' method description
#'
#' @param gplist ...
#' @param head ...
#' @param scales ...
#' @param highlight ...
#' @param ... Further arguments passed to method
#' @S3method xyplot gp.list
#' @method xyplot gp.list
xyplot.gp.list <- function(gplist, head = 30,
                           scales = list(x = list(rot = 90)),
                           highlight = NULL, ...){
  gpdata <- prepareGPPlot(gplist)

  if (!"Exact" %in% gpdata[, "method"]){
    subsetBy <- names(gplist)[1]
    warning(paste("Subsetting by", subsetBy))
  }

  # compute the probability of each graph
  subsetTotals <- tapply(gpdata[, c("gp")], gpdata[, "name"], sum)
  # sort the names by decreasing order of probability
  subsetNames <- names(sort(subsetTotals, dec = T))
  # sort the name factor by graph probability
  gpdata$name <- with(gpdata, reorder(name, gp))

  # reverse order of the levels
  gpdata$name <- factor(as.character(gpdata$name),
                        levels = rev(levels(gpdata$name)))

  if (!is.null(head)){
    subsetNames <- subsetNames[seq_len(head)]
  }

  numberOfTypes <- length(gplist)
  typesSeq <- seq_len(numberOfTypes)
  if (is.null(highlight)){
    gpdata$highlight <- factor(rep(T, nrow(gpdata)), levels = c(T, F))
  }
  else {
    gpdata$highlight <- factor(gpdata$name %in% highlight, levels = c(T, F))
  }

  typeSeq <- seq_len(numberOfTypes)
  mysuperpose <- trellis.par.get("superpose.symbol")
  mysuperpose$col <- rep(trellis.par.get("superpose.symbol")$col[typeSeq], 2)
  mysuperpose$pch <- rep(seq_len(numberOfTypes) + 2, 2)
  mysuperpose$alpha <- rep(c(1, 0.2), each = numberOfTypes)

  xyplot(
    gp ~ name,
    data               = gpdata,
    groups             = interaction(gpdata$method, highlight),
    type               = c("p", "h"),
    xlab               = "Graphs",
    ylab               = "Graph probability",
    scales             = scales,
    subset             = gpdata$name %in% subsetNames,
    key                = list(
                           text   = list(levels(gpdata$method)),
                           points = Rows(mysuperpose, typesSeq)
                         ),
    par.settings       = list(
      superpose.symbol = mysuperpose
    ),
    ...
  )
}

#' method name
#'
#' method description
#'
#' @param object ...
#' @param ... Further arguments passed to method
#' @S3method summary gp
#' @method summary gp
summary.gp <- function(object, ...){
  data.frame(
    graph = as.character(as.parental(names(object)), pretty = T),
    gp = as.numeric(object)
  )
}

#' Return the parental given by thresholding a edge probability matrix
#' at a given level. The inequality is >=. Note this may well be cyclic.
#'
#' @param ep A matrix of class ep
#' @param threshold The value at which to threshold the matrix
#' @return An object of class parental
#' @export
parentalFromEPThreshold <- function(ep, threshold){
  stopifnot("ep"              %in% class(ep),
            class(threshold)  ==   "numeric",
            length(threshold) ==   1)
  n <- dim(ep)[1]

  # threshold the edge probability table at threshold
  edgelist <- which(ep >= threshold, arr.ind = T)

  # remove self-loops
  edgelist <- edgelist[edgelist[, 1] != edgelist[, 2], , drop = F]

  # convert the resulting edgelist to parental graph list thing
  as.parental(edgelist, type = "edgelist", n)
}

#' Convert parental.list to CPDAGs
#'
#' A wrapper for \code{\link[parental]{as.cpdag}} that routes around the
#' issue that the parental.list is not of class \code{bn.list} etc.
#'
#' @param x A \code{parental.list}
#' @param verbose A logical. Should a progress bar be displayed?
#' @seealso \code{\link[parental]{as.cpdag}}, the function for which this
#'   function is a wrapper.
#' @return A \code{parental.list} of CPDAGs.
#' @export
parentalToCPDAG <- function(x, verbose = T){
  stopifnot(inherits(x, "parental.list"))

  # convert to bn, otherwise as.cpdag won't work.
  bnlist <- lapply(x, function(y){
    class(y) <- c("bn", "parental")
    y
  })
  class(bnlist) <- c("bn.list", "parental.list")

  if (verbose){
    progress <- create_progress_bar("text")
    progress$init(length(bnlist))
  }

  bnlistcp <- lapply(bnlist, function(y){
    if (verbose){
      progress$step()
    }
    as.cpdag(y)
  })

  # convert back to parental.list for ep(x, method = "tabulate")
  class(bnlistcp) <- "parental.list"
  bnlistcp
}

#' Prepare data for plotting in a ROC plot
#'
#' Converts MCMC samples, edge probability matrices etc to a data frame
#' ready for plotting as ROC curves.
#'
#' @param x The object to be prepared for an ROC plot
#' @param ... Further arguments passed to method
#' @seealso For actual graphs, see \code{\link{as.roc.parental}} and
#'   \code{\link{as.roc.parental.list}}. For edge probability matrices see
#'   \code{\link{as.roc.ep}} and \code{\link{as.roc.ep.list}}.
#' @return A data frame, with columns \code{estimate}, \code{tp}, and
#'   \code{fp}. The first contains the supplied label; the latter two
#'   contain the number of true and false positives respectively.
#' @export
as.roc <- function(x, ...){
  UseMethod("as.roc")
}

#' Prepare a parental for a ROC plot
#'
#' Compares two graphs, one of which is the true graph, and computes the
#' number of true and false positives.
#'
#' @param x A \code{parental} graph to compare to \code{true}
#' @param true A \code{bn}. The true graph.
#' @param label A label
#' @param ... Further arguments (unused)
#' @seealso For lists of graphs, see \code{\link{as.roc.parental.list}}. For
#'   edge probability matrices see \code{\link{as.roc.ep}} and
#'   \code{\link{as.roc.ep.list}}.
#' @return A 1-row data frame, with columns \code{estimate}, \code{tp}, and
#'   \code{fp}. The first contains the supplied label; the latter two
#'   contain the number of true and false positives respectively
#' @S3method as.roc parental
#' @method as.roc parental
as.roc.parental <- function(x, true, label, ...){
  tp <- pintersect(x, true, count = T)
  fp <- psetdiff(x, true, count = T)

  data.frame(estimate = label,
             tp   = tp,
             fp   = fp)
}

#' Prepare a parental list for a ROC plot
#'
#' Compares a list of graphs \code{x} to a true graph \code{true}, and
#' returns the number of true and false positives.
#'
#' @param x The object to be prepared for an ROC plot
#' @param true A \code{bn}. The true graph.
#' @param labels A vector of labels
#' @param verbose A logical. Should progress be displayed?
#' @param ... Further arguments (unused)
#' @seealso For individual graphs, see \code{\link{as.roc.parental}}. For
#'   edge probability matrices see \code{\link{as.roc.ep}} and
#'   \code{\link{as.roc.ep.list}}.
#' @return A data frame, with columns \code{estimate}, \code{tp}, and
#'   \code{fp}. The first contains the supplied label; the latter two
#'   contain the number of true and false positives respectively
#' @S3method as.roc parental.list
#' @method as.roc parental.list
as.roc.parental.list <- function(x, true, labels, verbose, ...){
  xSeq <- seq_along(x)
  out <- lapply(xSeq, function(i){
    as.roc(x[[i]], true, labels[[i]])
  })
  do.call("rbind", out)
}

#' Prepare an edge probability matrix for a ROC plot
#'
#' Compares an edge probability matrices \code{x} to a true graph
#' \code{true}, and returns the number of true and false positives.
#'
#' @param x A matrix of class \code{ep} of edge probabilities
#' @param true A \code{bn}. The true graph.
#' @param thresholds A numeric vector of thresholds.
#' @param label A label
#' @param verbose A logical. Should progress should be displayed?
#' @param ... Further arguments (unused)
#' @seealso For lists of edge probability matrices, see
#'   \code{\link{as.roc.ep.list}}. For graphs, \code{\link{as.roc.parental}}
#'   and \code{\link{as.roc.parental.list}}.
#' @return A data frame, with columns \code{estimate}, \code{tp}, and
#'   \code{fp}. The first contains the supplied label; the latter two
#'   contain the number of true and false positives respectively. Rows
#'   correspond to the supplied thresholds.
#' @S3method as.roc ep
#' @method as.roc ep
as.roc.ep <- function(x, true, thresholds, label, verbose, ...){
  rocdata <- data.frame(estimate = c(), tp = c(), fp = c())

  for (threshold in thresholds){
    xgraph <- parentalFromEPThreshold(x, threshold)
    newdata <- as.roc(xgraph, true, label = label)
    rocdata <- rbind(rocdata, newdata)
  }
  rocdata
}

#' Prepare an list of edge probability matrix for a ROC plot
#'
#' Compares a list of edge probability matrices \code{x} to a true graph
#' \code{true}, and returns the number of true and false positives.
#'
#' @param x A list of matrices of class \code{ep} of edge probabilities
#' @param true A \code{bn}. The true graph.
#' @param thresholds A numeric vector of thresholds.
#' @param labels A vector of labels
#' @param verbose A logical. Should progress should be displayed?
#' @param ... Further arguments (unused)
#' @seealso For individual edge probability matrices see
#'   \code{\link{as.roc.ep}}. For graphs, \code{\link{as.roc.parental}} and
#'   \code{\link{as.roc.parental.list}}.
#' @return A data frame, with columns \code{estimate}, \code{tp}, and
#'   \code{fp}. The first contains the supplied label; the latter two
#'   contain the number of true and false positives respectively. Rows
#'   correspond to the supplied thresholds, for each element of the supplied
#'   \code{ep.list}.
#' @S3method as.roc ep.list
#' @method as.roc ep.list
as.roc.ep.list <- function(x, true, thresholds, labels, verbose, ...){
  xSeq <- seq_along(x)
  out <- lapply(xSeq, function(i){
    as.roc(x[[i]], true, thresholds, labels[[i]])
  })
  do.call("rbind", out)
}

#' Plot an ROC curve
#' 
#' Print a ROC curve given a variety of estimation methods.
#' 
#'
#' An alternative is given by:
#' Werhli et al. Comparative evaluation of reverse engineering gene
#' regulatory networks with relevance networks, graphical gaussian models
#' and bayesian networks. Bioinformatics (2006) vol. 22 (20) pp. 2523
#'
#' @param true An object of class bn giving the true graph. This is
#'   converted to a CPDAG before comparision.
#' @param maps An object of class \code{bn.list}. An estimate of the true
#'   graph, given by the MAP estimate given by the MCMC. This is internally
#'   converted to CPDAG.
#' @param bnpmls An object of class \code{bnpostmcmc.list}. An object
#'   encapsulating a number of MCMC runs.
#' @param eps An object of class \code{ep.list}. 
#' @param thresholds An object of class numeric, giving the thresholds at
#'   which to plot the graphs given by MCMC and PC-MCMC estimates
#' @param use.cpdags A logical of length 1 indicating whether DAG should be 
#'   converted to CPDAGs before computing the ROC curves.
#' @param verbose A logical indicating whether to show progress bars etc.
#' @param ... Further arguments
#' @return A lattice object, containing the ROC plot
#' @export
rocplot <- function(true,
                    maps,
                    bnpmls,
                    eps,
                    thresholds = seq(from = 0, to = 1, by = 0.01),
                    use.cpdags = F,
                    verbose    = F,
                    ...){
  haveMAPs <- !missing(maps)
  haveBNPMLs <- !missing(bnpmls)
  haveEPs <- !missing(eps)
  
  stopifnot("bn" %in% class(true),
            !haveMAPs   || inherits(maps, "bn.list"),
            !haveBNPMLs || inherits(bnpmls, "bnpostmcmc.list"),
            !haveEPs    || inherits(eps, "ep.list"),
            inherits(thresholds, c("numeric", "integer")),
            inherits(use.cpdags, "logical"),
            inherits(verbose, "logical"))
  if (haveMAPs){
    mapsNames <- names(maps)
  }
  if (haveBNPMLs){
    bnpmlsNames <- names(bnpmls)
  }
  if (haveEPs){
    epsNames <- names(eps)
  }

  if (use.cpdags){
    if (verbose) cat("Converting true to cpdag\n")
    true <- as.cpdag(true)
    
    if (haveMAPs){
      if (verbose) cat("Converting MAPs to cpdag\n")
      maps <- lapply(maps, as.cpdag)
    }
    
    if (haveBNPMLs){
      if (verbose) cat("Converting BNPMLs to cpdag\n")
      bnpmls <- ep(x       = bnpmls,
                   method  = "tabulate",
                   FUN     = parentalToCPDAG,
                   verbose = T)
    }
  } else {
    if (haveBNPMLs){
      if (verbose) cat("Tabulating BNPMLs\n")
      bnpmls <- ep(x       = bnpmls,
                   method  = "tabulate",
                   verbose = T)
    }
  }

  rocdata <- data.frame(estimate = c(), tp = c(), fp = c())
  thresholds <- rev(thresholds)

  if (haveMAPs){
    if (verbose) cat("Computing ROC for Map\n")
    rocdataMap <- as.roc(x       = maps,
                         true    = true,
                         labels  = mapsNames,
                         verbose = verbose)
    rocdata <- rbind(rocdata, rocdataMap)
  }

  if (haveBNPMLs){
    if (verbose) cat("Computing ROC for BNPMLs\n")
    rocdataBNPMLs <- as.roc(x          = bnpmls,
                            true       = true,
                            thresholds = thresholds,
                            labels     = bnpmlsNames,
                            verbose    = verbose)
    rocdata <- rbind(rocdata, rocdataBNPMLs)
  }

  if (haveEPs){
    if (verbose) cat("Computing ROC for EPs\n")
    rocdataEP <- as.roc(x          = eps,
                        true       = true,
                        thresholds = thresholds,
                        labels     = epsNames,
                        verbose    = verbose)
    rocdata <- rbind(rocdata, rocdataEP)
  }

  rocdata[, "estimate"] <- factor(rocdata[, "estimate"])
  
  trimName <- function(x){
    lens <- sapply(x, nchar)
    substr(x, 1, lens - 2)
  }
  
  typeChar <- as.character(rocdata[, "estimate"])
  rocdata[, "type"] <- trimName(typeChar)
  ll <- c("Gibbs", "M-H", "REV", "Xie")
  rocdata[, "type"] <- factor(rocdata[, "type"], levels = ll)

  scalesAt <- c(0, seq_len(length(true) * (length(true) - 1)))

  xyplot(tp ~ fp,
         data         = rocdata,
         groups       = rocdata$type,
         xlab         = "False positives",
         ylab         = "True positives",
         type         = c("s", "s", "s", "p"),
         scales       = list(at = scalesAt),
         auto.key     = T,
         par.settings = simpleTheme(pch = 3:10),
         panel = panel.superpose,
         panel.groups = function(x, y, ...){
           panel.xyplot(x, y, groups = rocdata$estimate, ...)
         },
         ...)
}

#' method name
#'
#' method description
#'
#' @param exactgp ...
#' @param bnpostmcmclist ...
#' @param start ...
#' @param end ...
#' @param nbin ...
#' @param ... Further arguments passed to method
#' @export
cumtvd <- function(exactgp, bnpostmcmclist, start = 1, end,
                   nbin = floor((start-end)/100)){
  stopifnot(class(bnpostmcmclist) == "bnpostmcmc.list")
  if (missing(end)){
    end <- length(bnpostmcmclist[[1]]$samples)
  }

  numberOfNodes <- length(bnpostmcmclist[[1]]$samples[[1]])
  numberOfRuns <- length(bnpostmcmclist)

  # heuristic for whether to use pretty printing
  pretty <- seemsPretty(names(exactgp[1]))

  # get the graph probs
  # for each graph
  # for each MCMC sample
  # for each bin
  gpl <- gp(
    bnpostmcmclist,
    start = start,
    end = end,
    nbin = nbin,
    levels = names(exactgp),
    pretty = pretty
  )

  # want to convert the gp list to a matrix
  # with individual graphs in the columns
  # with the graphs ordered according to as.table in xyplot
  convertToMatrix <- function(gp){
    matrix(unlist(gp), ncol = numberOfGraphsDiscovered, byrow = T)
  }
  numberOfGraphsDiscovered <- length(gpl[[1]][[1]])
  gpmxl <- lapply(gpl, convertToMatrix)

  # convert to cumulative probabilities
  convertToCumulativeProbabilities <- function(gpmx, nbin){
    1/(seq_len(nbin)) * apply(gpmx, 2, cumsum)
  }
  cumgpmxl <- lapply(gpmxl, convertToCumulativeProbabilities, nbin)

  # convert the exact gp to a dataframe
  # with graphs across cols
  # and with the nbin number of rows
  exactgpmx <- matrix(exactgp, ncol = length(exactgp), nrow = nbin, byrow = T)

  # compute the tvd
  # for each MCMC sample
  # for each each bin
  cumtvds <- lapply(cumgpmxl, function(cumgpmx){
    rowSums(abs(cumgpmx - exactgpmx))
  })

  # do some nice naming
  names(cumtvds) <- paste("Run", seq_len(numberOfRuns))

  # return list of probs
  out <- list(
    cumtvds = cumtvds,
    lengthOfRuns = end - (start - 1),
    nbin = nbin,
    numberOfNodes = numberOfNodes,
    numberOfRuns = numberOfRuns
  )

  class(out) <- "cumtvd"
  out
}

#' method name
#'
#' method description
#'
#' @param cumtvd ...
#' @S3method xyplot cumtvd
#' @method xyplot cumtvd
xyplot.cumtvd <- function(cumtvd){
  # stack the cumtvd to a dataframe
  # with a column 'ind'
  # denoting which MCMC sample it orginated from
  data <- stack(cumtvd$cumtvd)

  # index along the bins
  data[[".index"]] <- seq_len(cumtvd$nbin)

  xyplot(
    values ~ .index,
    data,
    panel = panel.superpose,
    groups = data$ind,
    type = "l",
    main = paste(
      "Total variation distance through time, with",
      cumtvd$lengthOfRuns,
      "samples in",
      cumtvd$numberOfRuns,
      "runs"
    ),
    auto.key = T,
    ylab = "Total variation distance",
    xlab = "Samples",
    scales = list(x = list(draw = F))
  )
}

#' method name
#'
#' method description
#'
#' @param ... Further arguments passed to method
#' @export
cumtvd.list <- function(...){
  cumtvdlist <- list(...)
  class(cumtvdlist) <- c("cumtvd.list")
  cumtvdlist
}

#' method name
#'
#' method description
#'
#' @param cumtvdlist ...
#' @S3method xyplot cumtvd.list
#' @method xyplot cumtvd.list
xyplot.cumtvd.list <- function(cumtvdlist){
  cumtvdl <- lapply(cumtvdlist, function(cumtvd){
    # stack the cumtvd to a dataframe
    # with a column 'ind'
    # denoting which MCMC sample it orginated from
    stack(cumtvd$cumtvd)
  })
  names(cumtvdl) <- names(cumtvdlist)

  data <- do.call("make.groups", cumtvdl)

  # index along the bins
  data[[".index"]] <- seq_len(cumtvdlist[[1]]$nbin)

  xyplot(
    values ~ .index | which,
    data,
    groups = data$ind,
    type = "l",
    main = paste(
      "Total variation distance through time, with (for type 1)",
      cumtvdlist[[1]]$lengthOfRuns,
      "samples in",
      cumtvdlist[[1]]$numberOfRuns,
      "runs"
    ),
    auto.key = T,
    ylab = "Total variation distance",
    xlab = "Samples",
    scales = list(x = list(draw = F))
  )
}
