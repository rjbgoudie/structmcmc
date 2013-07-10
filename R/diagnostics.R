# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Edge probabilties matrix.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{epmx.bnpostmcmc.list}}, \code{\link{splom.epmx}},
#'   \code{\link{splom.bnpostmcmc.list}}
epmx <- function(x, ...){
  UseMethod("epmx")
}

#' Edge probabilities matrix.
#'
#' Computes the edge probabilities and return a matrix with these. The
#' format of the matrix is designed for the plotting function splom.epmx.
#'
#' For a problem with k nodes, the output will have k^2 columns and nbin
#' rows. Columns are in order 1->1, 1->2, 1->3, ...., 2->1, 2->2 etc
#'
#' @param x An object of class "bnpostmcmc.list"
#' @param method Either \code{"online"}, or \code{"offline"}. For online,
#'   the edge totals kept online will be used. This
#' @param nbin The number of equally-sized bins into which the samples are
#'          divided into before computing the edge probabilities of each.
#'          This applies only for \code{method = "offline"}.
#' @param start An integer of length 1, specifying the amount of burn-in.
#'          The samples start:end inclusive will be used.
#' @param end An integer of length 1, specifying the number of samples at the
#'          end to ignore. The samples start:end inclusive will be used.
#' @param verbose Should progress be shown? A logical.
#' @param ... Further arguments (unused)
#'
#' @return An object of class "epmx", a matrix of the form described above.
#' @S3method epmx bnpostmcmc.list
#' @method epmx bnpostmcmc.list
#' @seealso \code{\link{epmx}}, \code{\link{splom.epmx}},
#'   \code{\link{splom.bnpostmcmc.list}}
epmx.bnpostmcmc.list <- epmx.list <- function(x,
                                 method  = "online",
                                 nbin    = floor((end - start + 1)/100),
                                 start   = 1,
                                 end     = length(x[[1]]),
                                 verbose = F,
                                 ...){
  stopifnot(class(x) == "bnpostmcmc.list" || class(x) == "list",
            all(diff(sapply(x, length)) == 0) || class(x) == "list",
            method %in% c("online", "offline"),
            #validStartEnd(start, end, length(x[[1]])),
            is.wholenumber(nbin))
  haveEPL <- ifelse(class(x) == "list", T, F)

  if (haveEPL){
    numberOfNodes <- dim(x[[1]][[1]])[1]
    epl <- x
    epmxl <- lapply(epl, convertEPMatrixToColumns, numberOfNodes, verbose)
  }
  else {
    numberOfNodes <- length(x[[1]]$samples[[1]])
    if (method == "online"){
      samplers <- lapply(x, "[[", "sampler")
      epmxl <- lapply(samplers, epmx)
      nbin <- nrow(epmxl[[1]])
    } else {
      epmxl <- lapply(x, epmx, nbin, start, end, verbose, ...)
    }
  }

  # name the components
  names(epmxl) <- paste("Sample", seq_along(x))

  # return list of probs
  attr(epmxl, "lengthOfRuns") <-  end - (start - 1)
  attr(epmxl, "nbin") <- nbin
  attr(epmxl, "numberOfNodes") <- numberOfNodes
  attr(epmxl, "numberOfRuns") <- length(x)
  attr(epmxl, "type") <- "bnpostmcmc.list"

  class(epmxl) <- "epmx"
  epmxl
}


#' Edge probabilities matrix.
#'
#' Computes the edge probabilities and return a matrix with these. The
#' format of the matrix is designed for the plotting function splom.epmx.
#'
#' For a problem with k nodes, the output will have k^2 columns and nbin
#' rows. Columns are in order 1->1, 1->2, 1->3, ...., 2->1, 2->2 etc
#'
#' @param x An object of class "samplers"
#' @param verbose Should progress be shown? A logical.
#' @param ... Further arguments (unused)
#'
#' @return An object of class "epmx", a matrix of the form described above.
#' @S3method epmx samplers
#' @method epmx samplers
#' @seealso \code{\link{epmx}}, \code{\link{splom.epmx}},
#'   \code{\link{splom.bnpostmcmc.list}}
epmx.samplers <- function(x,
                          verbose = F,
                          ...){
  stopifnot(inherits(x, "samplers"))

  epmxl <- lapply(x, epmx)
  nbin <- nrow(epmxl[[1]])

  nSteps <- get("nSteps", envir = environment(x[[1]]))
  burnin <- get("nBurnin", envir = environment(x[[1]]))
  numberOfNodes <- get("numberOfNodes", envir = environment(x[[1]]))

  # name the components
  names(epmxl) <- paste("Sample", seq_along(x))

  # return list of probs
  attr(epmxl, "lengthOfRuns") <- nSteps
  attr(epmxl, "nbin") <- nbin
  attr(epmxl, "numberOfNodes") <- numberOfNodes
  attr(epmxl, "numberOfRuns") <- length(x)
  attr(epmxl, "type") <- "bnpostmcmc.list"

  class(epmxl) <- "epmx"
  epmxl
}

#' Edge probabilities matrix from a sampler.
#'
#' Computes the edge probabilities and return a matrix with these. The
#' format of the matrix is designed for the plotting function splom.epmx.
#'
#' For a problem with k nodes, the output will have k^2 columns and nbin
#' rows. Columns are in order 1->1, 1->2, 1->3, ...., 2->1, 2->2 etc
#'
#' @param x An MCMC sampler, of class "sampler".
#' @param verbose Should progress be shown? A logical.
#' @param ... Further arguments (unused)
#'
#' @return An object of class "epmx", a matrix of the form described above.
#' @S3method epmx sampler
#' @method epmx sampler
#' @seealso \code{\link{bnpostmcmc.list}}, \code{\link{epmx}},
#'   \code{\link{splom.epmx}}, \code{\link{splom.bnpostmcmc.list}}
epmx.sampler <- function(x, verbose = F, ...){
  stopifnot(inherits(x, "sampler"))

  etbins_exists <- exists("etbins", envir = environment(x))

  if (etbins_exists){
    epmx <- get("etbins", envir = environment(x))
    epmx <- epmx[!is.na(rowSums(epmx)), , drop = F]

    # compute edge probabilities from totals
    nSteps <- get("nSteps", envir = environment(x))
    nBurnin <- get("nBurnin", envir = environment(x))
    etBinsSize <- get("etBinsSize", envir = environment(x))
    epmx <- epmx/etBinsSize

    # adjust for last row being incomplete
    leftover <- (nSteps - nBurnin) %% etBinsSize
    if (leftover > 0){
      epmx[nrow(epmx), ] <- (epmx[nrow(epmx), ] * etBinsSize)/leftover
    }
    epmx
  } else {
    stop("et bins missing")
  }
}

#' Convert edge prob matrix to a column matrix
#'
#' Want to convert the ep list to a matrix with individual edges in the
#' columns with the edges ordered according to as.table in splom
#'
#' ie 1->1, 1->2, 1->3, ... , 2->1, 2->2, ...
#' @param epl A matrix of edge probability matrices, each corresponding to
#'   a separate bin
#' @param numberOfNodes The number of nodes in the matrix
#' @param verbose Should progress text be output?
#' @export
convertEPMatrixToColumns <- function(epl, numberOfNodes, verbose){
  if (verbose){
    cat("Converting to Matrix\n")
  }
  epl <- lapply(epl, t)
  if (verbose){
    cat("matrix line\n")
  }
  matrix(unlist(epl, use.names = F), ncol = numberOfNodes^2, byrow = T)
}

#' Edge probabilities matrix.
#'
#' Computes the edge probabilities and return a matrix with these. The
#' format of the matrix is designed for the plotting function splom.epmx.
#'
#' For a problem with k nodes, the output will have k^2 columns and nbin
#' rows. Columns are in order 1->1, 1->2, 1->3, ...., 2->1, 2->2 etc
#'
#' @param x An object of class "bnpostmcmc.list"
#' @param nbin The number of equally-sized bins into which the samples are
#'          divided into before computing the edge probabilities of each
#' @param start An integer of length 1, specifying the amount of burn-in.
#'          The samples start:end inclusive will be used.
#' @param end An integer of length 1, specifying the number of samples at the
#'          end to ignore. The samples start:end inclusive will be used.
#' @param verbose Should progress be shown? A logical.
#' @param ... Further arguments (unused)
#'
#' @return An object of class "epmx", a matrix of the form described above.
#' @S3method epmx bnpostmcmc
#' @method epmx bnpostmcmc
#' @seealso \code{\link{bnpostmcmc.list}}, \code{\link{epmx}},
#'   \code{\link{splom.epmx}}, \code{\link{splom.bnpostmcmc.list}}
epmx.bnpostmcmc <- function(x,
                            nbin    = floor((end - start + 1)/100),
                            start   = 1,
                            end     = length(x[[1]]),
                            verbose = F,
                            ...){
  stopifnot(class(x) == "bnpostmcmc",
            is.wholenumber(nbin))

  numberOfNodes <- length(x$samples[[1]])
  eps <- ep(x, start = start, end = end, nbin = nbin, method = "flatten")

  convertEPMatrixToColumns(eps, numberOfNodes, verbose)
}

#' Cumulative mean.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{cummean.matrix}}
cummean <- function(x, ...){
  UseMethod("cummean")
}

#' Cumulative mean.
#'
#' Compute the cumulative means of the columns of a matrix x.
#' ie each column is treated separately
#'
#' @param x A matrix
#' @param ... Further arguments (unused)
#'
#' @return A matrix of the same dimension as x, with the cumulative values.
#' @S3method cummean matrix
#' @method cummean matrix
#' @seealso \code{\link{cummean}}
cummean.matrix <- function(x, ...){
  result <- 1/(seq_len(nrow(x))) * apply(x, 2, cumsum)
  dim(result) <- dim(x)
  result
}

#' Cumulative edge probabilities.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{cumep.bnpostmcmc.list}}
cumep <- function(x, ...){
  UseMethod("cumep")
}

#' Cumulative edge probabilities.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @S3method cumep bnpostmcmc.list
#' @S3method cumep list
#' @method cumep bnpostmcmc.list
#' @method cumep list
#' @seealso \code{\link{cumep}}
cumep.bnpostmcmc.list <- cumep.list <- function(x, ...){

  res <- epmx(x, ...)
  at <- attributes(res)
  res <- lapply(res, cummean)
  attributes(res) <- at
  attr(res, "function") <- "cum"
  res
}

#' Cumulative edge probabilities.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @S3method cumep samplers
#' @method cumep samplers
#' @seealso \code{\link{cumep}}
cumep.samplers <- function(x, ...){
  res <- epmx(x, ...)
  at <- attributes(res)
  res <- lapply(res, cummean)
  attributes(res) <- at
  attr(res, "function") <- "cum"
  res
}

#' Moving window mean.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{mwmean.matrix}}
mwmean <- function(x, ...){
  UseMethod("mwmean")
}

#' Moving window mean.
#'
#' Compute the moving window means of the columns of a matrix x.
#' ie each column is treated separately
#'
#' @param x A matrix
#' @param window ...
#' @param ... Further arguments (unused)
#'
#' @return A matrix of the same dimension as x, with the moving window values.
#' @S3method mwmean matrix
#' @method mwmean matrix
#' @seealso \code{\link{mwmean}}
mwmean.matrix <- function(x, window, ...){
  if (require("zoo")){
    unname(as.matrix(rollapply(as.zoo(x), window, mean)))
  }
}

#' Moving window edge probiilities.
#'
#' method description
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{mwep.bnpostmcmc.list}}
mwep <- function(x, ...){
  UseMethod("mwep")
}

#' Moving window edge probabilities.
#'
#' method description
#'
#' @param x ...
#' @param window ...
#' @param ... Further arguments passed to method
#' @S3method mwep bnpostmcmc.list
#' @S3method mwep list
#' @method mwep bnpostmcmc.list
#' @method mwep list
#' @seealso \code{\link{mwep}}
mwep.bnpostmcmc.list <- mwep.list <- function(x, window = 10, ...){
  res <- epmx(x, ...)
  at <- attributes(res)
  res <- lapply(res, mwmean, window)
  attributes(res) <- at
  attr(res, "function") <- "mw"
  res
}

#' Moving window edge probabilities.
#'
#' method description
#'
#' @param x ...
#' @param window ...
#' @param ... Further arguments passed to method
#' @S3method mwep samplers
#' @method mwep samplers
#' @seealso \code{\link{mwep}}
mwep.samplers <- function(x, window = 10, ...){
  res <- epmx(x, ...)
  at <- attributes(res)
  res <- lapply(res, mwmean, window)
  attributes(res) <- at
  attr(res, "function") <- "mw"
  res
}

#' Plot of cumulative edge probabilities.
#'
#' Returns a splom of the cumulative edge probabilities through time
#' for bnpostmcmc.list and bvspostmcmc.list
#'
#' @param x An \code{epmx} object.
#' @param subset A numeric vector that specifies the nodes between which
#'   the plot will be drawn for. The default value \code{NULL} displays
#'   all pairs.
#'
#' @return A scatter plot matrix.
#' @S3method splom epmx
#' @method splom epmx
#' @seealso \code{\link{epmx}}
splom.epmx <- function(x, subset = NULL){
  if (!is.null(subset)){
    stopifnot(inherits(subset, c("numeric", "integer")))
    subset <- list(from = subset, to = subset)
  }
  epmxPlotInternal(x, subset, plottype = "splom")
}

#' Plot of cumulative edge probabilities.
#'
#' Returns a xyplot of the cumulative edge probabilities through time
#' for bnpostmcmc.list and bvspostmcmc.list
#'
#' @param x An \code{epmx} object.
#' @param subset A list of length of two, the first component of which
#'   determines the heads of the edges that are displayed, and the second
#'   determines the tails of the edges that are displayed. The default
#'   value \code{NULL} displays all pairs.
#'
#' @return An \code{xyplot}.
#' @S3method xyplot epmx
#' @method xyplot epmx
#' @seealso \code{\link{epmx}}
xyplot.epmx <- function(x, subset = NULL){
  epmxPlotInternal(x, subset, plottype = "xyplot")
}

#' Convert epmx to a data.frame
#'
#' Converts an epmx to a data.frame, in the appropriate form for
#' \code{\link{epmxPlotInternal}}.
#'
#' @param x An \code{epmx} object.
#' @param subset A list of length of two, the first component of which
#'   determines the heads of the edges that are displayed, and the second
#'   determines the tails of the edges that are displayed. The default
#'   value \code{NULL} displays all pairs.
#' @param plottype Either \code{"xyplot"} or \code{"splom"}.
#'
#' @return A \code{data.frame}, with a column \code{which} that denotes
#'   the sampler, a \code{.index} column indicating the original row.
#' @S3method as.data.frame epmx
#' @method as.data.frame epmx
#' @seealso \code{\link{epmx}}
as.data.frame.epmx <- function(x, subset, plottype = "xyplot"){
  numberOfNodes <- attr(x, "numberOfNodes")

  # reorder the columns to be appropriate for as.table
  # order is a matrix that converts colwise ordering to rowwise ordering
  dyadSeq <- seq_len(numberOfNodes^2)
  order <- matrix(dyadSeq, numberOfNodes, numberOfNodes, byrow = T)

  # Name the columns
  allDyads <- which(order > 0, arr.ind = T)
  newColnames <- apply(allDyads, 1, paste, collapse = "->")

  reorderAndAddColnames <- function(x){
    colnames(x) <- newColnames
    as.data.frame(x[, as.vector(order)])
  }

  x <- lapply(x, reorderAndAddColnames)

  # subsetting
  if (is.null(subset)){
    subset <- list(seq_len(numberOfNodes), seq_len(numberOfNodes))
  }

  stopifnot(inherits(subset, "list"),
            length(subset) == 2,
            all(sapply(subset[[1]], is.wholenumber)),
            all(sapply(subset[[2]], is.wholenumber)),
            all(findInterval(subset[[1]], c(1, numberOfNodes), r = T) == 1),
            all(findInterval(subset[[2]], c(1, numberOfNodes), r = T) == 1))

  # the subset[[2]] and subset[[1]] are transposed here because the
  # columns have already been reordered
  subsetm <- matrix(F, numberOfNodes, numberOfNodes)
  subsetm[subset[[2]], subset[[1]]] <- T
  diag(subsetm) <- F
  x <- lapply(x, "[", , as.vector(subsetm), drop = F)

  data <- do.call("make.groups", x)

  # add column denoting which bin the data is from
  s <- sapply(unlist(sapply(x, nrow)), seq)
  data <- cbind(data, `.index` = as.vector(unlist(s)))
  data
}

#' (Internal) Plot of cumulative edge probabilities.
#'
#' Returns a xyplot/splom of the cumulative edge probabilities through time
#' for bnpostmcmc.list and bvspostmcmc.list
#'
#' @param x An \code{epmx} object.
#' @param subset A list of length of two, the first component of which
#'   determines the heads of the edges that are displayed, and the second
#'   determines the tails of the edges that are displayed. The default
#'   value \code{NULL} displays all pairs.
#' @param plottype Either \code{"xyplot"} or \code{"splom"}.
#'
#' @return An appropriate plot.
#' @seealso \code{\link{epmx}}
epmxPlotInternal <- function(x, subset, plottype = "xyplot"){
  stopifnot(class(x) == "epmx")

  lengthOfRuns <- attr(x, "lengthOfRuns")
  numberOfRuns <- attr(x, "numberOfRuns")
  type <- attr(x, "type")
  isbvsresponse <- isTRUE(type == "bvspostmcmc.list")
  nbin <- attr(x, "nbin")
  fun <- attr(x, "function")

  data <- as.data.frame(x, subset = subset, plottype = plottype)

  # select all the columns apart from 'which'
  indexAndWhich <- (length(names(data)) - 1):length(names(data))
  variables <- names(data)[-indexAndWhich]
  variablesasname <- lapply(variables, as.name)
  toparse <- paste(paste(variablesasname, collapse = "+"), "~ .index")
  form <- eval(parse(text = toparse))

  if (plottype == "splom"){
    numberOfNodesToPlot <- length(subset[[1]])
  }
  mainintro <- switch(fun,
    cum = "Edge probabilities through time, with",
    mw = "Moving windows probs, with"
  )

  main <- paste(
    mainintro,
    lengthOfRuns,
    "samples in",
    numberOfRuns, "runs"
  )

  if (plottype == "xyplot"){
    xyplot(
      form,
      data = data,
      outer = T,
      groups = which,
      as.table = T,
      default.scales = list(relation = "same"),
      type = "l",
      ylab = "Probability",
      xlab = "Samples",
      ylim = c(-0.025, 1.025),
      scales = list(
        y = list(tick.number = 3, at = seq(from = 0, to = 1, by = 0.1)),
        x = list(draw = F)
      ),
      main = main,
      panel = function(x, y, varnames, ...){
        for (i in seq(from = 0, to = 1, by = 0.1)){
          panel.abline(h = i, col.line = rgb(0.9, 0.9, 0.9))
        }
        panel.xyplot(x, y, ...)
      },
      par.settings = list(
        axis.line = list(col = rgb(0.9, 0.9, 0.9))
      ),
      between = list(x = 0.5, y = 0.5)
    )
  } else {
    if (isbvsresponse){
      xyplot(
        form,
        data = data,
        outer = T,
        groups = which,
        as.table = T,
        default.scales = list(relation = "same"),
        type = "l",
        ylab = "Probability",
        xlab = "Samples",
        ylim = c(-0.025, 1.025),
        scales = list(
          y = list(tick.number = 3, at = seq(from = 0, to = 1, by = 0.1)),
          x = list(draw = F)
        ),
        main = main,
        panel = function(x, y, ...){
          for (i in seq(from = 0, to = 1, by = 0.1)){
            panel.abline(h = i, col.line = rgb(0.9, 0.9, 0.9))
          }
          panel.xyplot(x, y, ...)
        },
        par.settings = list(
          axis.line = list(col = rgb(0.9, 0.9, 0.9))
        ),
        between = list(x = 0.5, y = 0.5)
      )
    }
    else {
      xyplot(
        form,
        data = data,
        outer = T,
        groups = which,
        layout = c(numberOfNodesToPlot, numberOfNodesToPlot),
        as.table = T,
        default.scales = list(relation = "same"),
        strip = F,
        varnames = subset[[1]],
        type = "l",
        ylab = "Probability",
        xlab = "Samples",
        ylim = c(-0.025, 1.025),
        scales = list(
          y = list(tick.number = 3, at = seq(from = 0, to = 1, by = 0.1)),
          x = list(draw = F)
        ),
        main = main,
        panel = function(x, y, varnames, ...){
          if (current.row() == current.column()){
            diag.panel.splom(varname = varnames[current.row()],
                             limits  = c(0, 1),
                             draw    = F)
          }
          else {
            for (i in seq(from = 0, to = 1, by = 0.1)){
              panel.abline(h = i, col.line = rgb(0.9, 0.9, 0.9))
            }
          panel.xyplot(x, y, ...)
          }
        },
        par.settings = list(
          axis.line = list(col = rgb(0.9, 0.9, 0.9))
        ),
        between = list(x = 0.5, y = 0.5)
      )
    }
  }
}

#' Scatterplot matrix of edge probabilities between runs.
#'
#' method description
#'
#' @param x ...
#' @param start ...
#' @param end ...
#' @param plot ...
#' @S3method splom bnpostmcmc.list
#' @S3method splom list
#' @method splom bnpostmcmc.list
#' @method splom list
#' @seealso \code{\link{epmx}}, \code{\link{splom.epmx}}, \code{\link{splom}}
splom.bnpostmcmc.list <- splom.list <- function(x, start = 1, end, plot = F){
  stopifnot(class(x) == "bnpostmcmc.list" ||
                class(x) == "list",
            length(x) > 1)
  haveEPL <- ifelse(class(x) == "list", T, F)
  if (haveEPL){
    numberOfNodes <- dim(x[[1]])[1]
    numberOfRuns <- length(x)
    epl <- x
  } else {
    if (missing(end)){
      end <- length(x[[1]]$samples)
    }
    numberOfNodes <- length(x[[1]]$samples[[1]])
    numberOfRuns <- length(x)
    # get a matrix of edge probabilities
    # for each MCMC sample
    epl <- ep(x, start, end)
  }

  # remove the diagonal elements
  # which are not of interest
  epl <- lapply(epl, notdiag)

  # set up names
  names(epl) <- paste("Run", seq_len(numberOfRuns))

  # pre-allocate data matrix
  data <- matrix(ncol = length(epl), nrow = numberOfNodes ^ 2)
  # convert epl to a dataframe
  data <- data.frame(epl)

  splom(
    data,
    main = paste("Final edge probabilities from ", numberOfRuns, "runs"),
    par.settings = list(
      axis.line = list(col = rgb(0.9, 0.9, 0.9))
    ),
    between = list(x = 0.5, y = 0.5),
    xlim = c(-0.1, 1.1),
    ylim = c(-0.1, 1.1),
    pscales = lapply(
      seq_len(numberOfRuns),
      function(i){
        list(limits = c(-0.1, 1.1), at = c(0, 0.5, 1))
      }
    ),
    scales = list(relation = "same"),
    superpanel = function(panel, upper.panel, lower.panel, ...){
      panel.pairs(
        panel = lattice.getOption("panel.splom"),
        upper.panel = function(...){
          panel.refline(a = 0, b = 1, col.line = "lightblue", lwd = 1)
          panel.refline(h = 0)
          panel.refline(v = 0)
          panel.refline(h = 1)
          panel.refline(v = 1)
          panel.xyplot(...)
        },
        lower.panel = function(...){},
        as.matrix = T,
        ...
      )
    }
  )
}

# # @S3method dotplot gp
# xyplot.gp <- function(samples, tabulated = NULL, rank = NULL, ...){
#   if (is.null(tabulated)){
#     tabulated <- tabulate.samples(samples, prettyPrint = T)
#   }
#
#   dat <- data.frame(probability = as.numeric(tabulated/length(samples)))
#
#   dat$name <- names(tabulated)
#   dat$name <- with(dat, reorder(factor(names(tabulated)), probability))
#
#   if (!is.null(rank)){
#     dat <- dat[seq_len(rank), ]
#   }
#
#   dotplot(
#     name ~ probability,
#     data = dat,
#     xlab = "Graph probability",
#     type = c("p", "h"),
#     scales = list(y = list(cex = 0.1))
#   )
# }

#' Plot tape.
#'
#' method description
#'
#' @param sampler ...
#' @param burnIn ...
#' @export
plotTape <- function(sampler, burnIn = 10000){
  # get the tape from the sampler
  tape <- sampler(returnTape = T)
  nSteps <- nrow(tape)

  # chop off burnIn, after sanity checking
  if (nSteps <= burnIn){
    stop("burnIn must be shorter than the run length.")
  }
  tape <- tape[seq(from = burnIn, to = nSteps), ]

  # convert to the correct form for lattice
  tape <- as.data.frame(tape)
  tape <- cbind(tape, ind = seq_len(nrow(tape)))

  xyplot(
    logAcceptanceProb ~ ind | equal.count(ind, 9, overlap = 0.1),
    data = tape,
    groups = interaction(tape$movetype, tape$logAcceptanceProb > 20),
    xlab = "Steps",
    ylab = "Log acceptance probability",
    main = paste("Log acceptance probabilities from", nSteps, "samples"),
    pch = c(".", ".", ".", "."),
    cex = c(1, 1, 2, 2),
    strip = F,
    layout = c(1, 9),
    par.settings = list(
      axis.line = list(col = rgb(0.9, 0.9, 0.9))
    ),
    between = list(x = 0.5, y = 0.5),
    scales = list(x = list(relation = "sliced", axs = "i", draw = F))
  )
}

#' Gelman and Rubin's convergence diagnostic
#'
#' The 'potential scale reduction factor' is calculated for each variable in
#' \code{x}, together with upper and lower confidence limits.
#' Computation is performed using \code{\link[coda]{gelman.diag}}.
#'
#' @param x ...
#' @param ... Further arguments passed to method
#' @export
#' @seealso \code{\link{gelman.samplers}}. Computation performed using
#'   \code{\link[coda]{gelman.diag}}.
gelman <- function(x, ...){
  UseMethod("gelman")
}

#' Gelman and Rubin's convergence diagnostic
#'
#' The 'potential scale reduction factor' is calculated for each variable in
#' \code{x}, together with upper and lower confidence limits.
#' Computation is performed using \code{\link[coda]{gelman.diag}}.
#'
#' @param x An MCMC sampler
#' @param names A character vector of statistics collected during the
#'   MCMC run, for which the 'potential scale reduction factor' should be
#'   computed.
#' @param transform A list of functions (or a single function), which will
#'   be applied to the statistics. This can be used to improve the
#'   normality of the statistics, which is a requirement for Gelman-Rubin's
#'   statistic to be meaningful. The default \code{NULL} value applies no
#'   transform to any statistic (equivalent to \code{\link{identity}})
#' @param ... Further arguments, currently unused
#'
#' @return An object of class \code{gelman.diag}
#' @S3method gelman samplers
#' @method gelman samplers
#' @seealso Computation performed using \code{\link[coda]{gelman.diag}}.
gelman.samplers <- function(x, names = "nEdges", transform = NULL, ...){
  stopifnot(inherits(x, "samplers"),
            inherits(names, "character"))
  if (!is.null(transform)){
    if (length(names) == 1){
      stopifnot(inherits(transform, "function"))
      transform <- list(transform)
    }
    stopifnot(length(transform) == length(names))
  }

  if (require(coda)){
    ml <- lapply(x, function(sampler){
      ml <- statistics(sampler, names)
      if (!is.null(transform)){
        for (i in seq_along(names)){
          ml[[i]] <- transform[[i]](ml[[i]])
        }
      }
      ml <- data.frame(ml)
      mcmc(ml)
    })
    ml <- as.mcmc.list(ml)
    gelman.diag(ml)$psrf[1, ]
  } else {
    stop("Package 'coda' required to compute Gelman-Rubin statistic")
  }
}
