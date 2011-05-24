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
#' @seealso \code{\link{epmx.bnpostmcmc.list}}, \code{\link{xyplot.epmx}},
#'   \code{\link{splom.bnpostmcmc.list}}
epmx <- function(x, ...){
  UseMethod("epmx")
}

#' Edge probabilities matrix.
#' 
#' Computes the edge probabilities and return a matrix with these. The
#' format of the matrix is designed for the plotting function xyplot.epmx.
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
#' @seealso \code{\link{epmx}}, \code{\link{xyplot.epmx}},
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
#' format of the matrix is designed for the plotting function xyplot.epmx.
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
#' @seealso \code{\link{epmx}}, \code{\link{xyplot.epmx}},
#'   \code{\link{splom.bnpostmcmc.list}}
epmx.samplers <- function(x,
                          verbose = F,
                          ...){
  stopifnot(inherits(x, "samplers"))

  epmxl <- lapply(x, epmx)
  nbin <- nrow(epmxl[[1]])

  nSteps <- get("nSteps", envir = environment(x[[1]]))
  burnin <- get("burnin", envir = environment(x[[1]]))
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
#' format of the matrix is designed for the plotting function xyplot.epmx.
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
#'   \code{\link{xyplot.epmx}}, \code{\link{splom.bnpostmcmc.list}}
epmx.sampler <- function(x, verbose = F, ...){
  stopifnot(inherits(x, "sampler"))

  etbins_exists <- exists("etbins", envir = environment(x))

  if (etbins_exists){
    epmx <- get("etbins", envir = environment(x))
    epmx <- epmx[!is.na(rowSums(epmx)), , drop = F]

    # compute edge probabilities from totals
    nSteps <- get("nSteps", envir = environment(x))
    etBinsSize <- get("etBinsSize", envir = environment(x))
    epmx <- epmx/etBinsSize

    # adjust for last row being incomplete
    leftover <- nSteps %% etBinsSize
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
#' columns with the edges ordered according to as.table in xyplot 
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
#' format of the matrix is designed for the plotting function xyplot.epmx.
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
#'   \code{\link{xyplot.epmx}}, \code{\link{splom.bnpostmcmc.list}}
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
  1/(seq_len(nrow(x))) * apply(x, 2, cumsum)
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
  if (require(zoo)){
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
#' Returns a xyplot of the cumulative edge probabilities through time
#' for bnpostmcmc.list and bvspostmcmc.list
#'
#' @param x ...
#'
#' @return ...
#' @S3method xyplot epmx
#' @method xyplot epmx
#' @seealso \code{\link{epmx}}
xyplot.epmx <- function(x){

  stopifnot(class(x) == "epmx")
  lengthOfRuns <- attr(x, "lengthOfRuns")
  numberOfRuns <- attr(x, "numberOfRuns")
  numberOfNodes <- attr(x, "numberOfNodes")
  type <- attr(x, "type")
  nbin <- attr(x, "nbin")
  fun <- attr(x, "function")
  isbvsresponse <- isTRUE(type == "bvspostmcmc.list")

  # pre-allocate the data matrix
  if (isbvsresponse){
    data <- matrix(nrow = numberOfRuns * nbin, ncol = numberOfNodes)
  }
  else {
    data <- matrix(nrow = numberOfRuns * nbin, ncol = numberOfNodes ^ 2)
  }

  # convert each epmx x to a dataframe
  xDF <- lapply(x, function(x){
    out <- as.data.frame(x)
    # would be good to name these
    #colnames(out) <-
  })

  # stack with a column called 'which'
  # that shows which sample it came from
  data <- do.call("make.groups", xDF)

  # add column denoting which bin the data is from
  data[[".index"]] <- as.vector(sapply(sapply(xDF, nrow), seq))

  # select all the columns apart from 'which'
  indexAndWhich <- (length(names(data)) - 1):length(names(data))
  variables <- names(data)[-indexAndWhich]
  variablesasname <- lapply(variables, as.name)
  toparse <- paste(paste(variablesasname, collapse = "+"), "~ .index")
  form <- eval(parse(text = toparse))

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
      layout = c(numberOfNodes, numberOfNodes),
      as.table = T,
      default.scales = list(relation = "same"),
      strip = F,
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
        if (current.row() == current.column()){
          diag.panel.splom(varname = current.row(),
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
#' @seealso \code{\link{epmx}}, \code{\link{xyplot.epmx}}, \code{\link{splom}}
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
        list(limits = c(-0.1, 1.1), at = c(5))
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
