# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Draw samples from a MCMC sampler.
#'
#' Draws samples from an MCMC sampler. The length of the run can be specified
#' either by the number of samples to be drawn, or the length of time that
#' the sampler runs.
#'
#' @param sampler A sampler object
#' @param n The number of samples to draw. Set this to \code{FALSE} if
#'   using the \code{time} argument.
#' @param time The number of seconds to spend drawing samples. Set this to
#'   \code{FALSE} if using the \code{n} argument.
#' @param burnin The number of samples to discard from the beginning of
#'   the sample.
#' @param thin The frequency with which samples should be kept. eg for
#'   \code{thin = 3}, every third sample will be kept.
#' @param verbose A logical. Should a progress bar be displayed?
#' @export
#' @seealso \code{\link{drawSamplesByTime}},
#'   \code{\link{drawSamplesByStepCount}}. \code{\link{BNSampler}},
#'   \code{\link{BNGibbsSampler}}, \code{\link{BNSamplerMJ}}
#' @examples
#' set.seed(310)
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#'
#' prior <- function(net) 1
#' initial <- bn(c(), c())
#'
#' sampler <- BNSampler(dat, initial, prior)
#' samples <- draw(sampler, n = 5)
draw <- function(sampler, n = F, time = F, burnin = 0, thin = 1, verbose = T){
  stopifnot(identical(n, F) || inherits(n, c("integer", "numeric")),
            identical(n, F) || is.wholenumber(n),
            identical(F, time) || inherits(time, c("integer", "numeric")),
            inherits(burnin, c("integer", "numeric")),
            is.wholenumber(burnin),
            inherits(verbose, "logical"),
            inherits(thin, "logical") ||
              (inherits(thin, c("integer", "numeric")) &&
               is.wholenumber(thin)))

  nSteps <- get("nSteps", envir = environment(sampler))

  if (nSteps < burnin){
    burninleft <- burnin - nSteps
    if (inherits(time, c("integer", "numeric")) && thin > 1){
      burninleft <- burninleft + thin - burninleft %% thin
    }
  } else {
    burninleft <- 0
  }

  runLengthByTime <- identical(n, F)
  if (runLengthByTime){
    drawSamplesByTime(sampler    = sampler,
                      time       = time,
                      burnin     = burnin,
                      burninleft = burninleft,
                      thin       = thin,
                      verbose    = verbose)
  } else {
    lengthout <- max(0, floor((n - burninleft)/thin))
    drawSamplesByStepCount(sampler    = sampler,
                           n          = n,
                           burnin     = burnin,
                           burninleft = burninleft,
                           lengthout  = lengthout,
                           thin       = thin,
                           verbose    = verbose)
  }
}


#' Draw samples from a MCMC sampler, by time.
#'
#' Draws samples from an MCMC sampler for a specified number of seconds.
#'
#' @param sampler A sampler object
#' @param time The number of seconds to spend drawing samples.
#' @param burnin The number of samples to discard from the beginning of
#'   the sample (includes previous runs)
#' @param burninleft The number of samples still to be discarded.
#' @param thin The frequency with which samples should be kept. eg for
#'   \code{thin = 3}, every third sample will be kept.
#' @param verbose A logical. Should a progress bar be displayed?
#' @param samplesIncrement The size by which the samples should be
#'   incremented.
#' @export
#' @seealso \code{\link{draw}}, \code{\link{drawSamplesByStepCount}}
#' @examples
#' set.seed(310)
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#'
#' prior <- function(net) 1
#' initial <- bn(c(), c())
#'
#' sampler <- BNSampler(dat, initial, prior)
#' samples <- draw(sampler, time = 15)
drawSamplesByTime <- function(sampler,
                              time,
                              burnin,
                              burninleft,
                              thin,
                              verbose,
                              samplesIncrement = 10000){
  # Note: the initial graph is NOT returned at the moment
  samples <- vector("list", samplesIncrement)
  if (verbose){
    cat("Drawing samples, for about", time, "seconds\n")
    flush.console()
    pb <- txtProgressBar(max = time, style = 3)
    setTxtProgressBar(pb, 0)
  }

  i <- 1
  elapsed <- 1
  start <- proc.time()
  while (elapsed < time + 1){

    out <- sampler(i, burnin = burnin)

    if (i %% thin == 0 && i > burninleft){
      class(out) <- c("bn", "parental")
      samples[[(i-burninleft)/thin]] <- out
    }

    if (i %% samplesIncrement == 0){
      oldLength <- length(samples)
      newLength <- oldLength + samplesIncrement
      newSamples <- vector("list", newLength)
      newSamples[1:oldLength] <- samples
      samples <- newSamples
    }

    elapsed <- (proc.time() - start)[[3]]
    if (verbose){
      setTxtProgressBar(pb, value = elapsed)
    }
    i <- i + 1
  }

  if (verbose){
    close(pb)
    cat("Drew", i, "samples in", elapsed, "seconds")
  }

  samples <- samples[seq_len(max(0,i-burninleft - 1))]
  class(samples) <- c("mcmcbn", "bn.list", "parental.list")
  samples
}

#' Draw samples from a MCMC sampler, by step count.
#'
#' Draws a specific number of samples from an MCMC sampler.
#'
#' @param sampler A sampler object
#' @param n The number of samples to draw.
#' @param burnin The number of samples to discard from the beginning of
#'   the sample (includes previous runs)
#' @param burninleft The number of samples still to be discarded.
#' @param lengthout The length of the final samples object
#' @param thin The frequency with which samples should be kept. eg for
#'   \code{thin = 3}, every third sample will be kept.
#' @param verbose A logical. Should a progress bar be displayed?
#' @export
#' @seealso \code{\link{draw}}, \code{\link{drawSamplesByTime}}
#' @examples
#' set.seed(310)
#' x1 <- factor(c(1, 1, 0, 1))
#' x2 <- factor(c(0, 1, 0, 1))
#' dat <- data.frame(x1 = x1, x2 = x2)
#'
#' prior <- function(net) 1
#' initial <- bn(c(), c())
#'
#' sampler <- BNSampler(dat, initial, prior)
#' samples <- draw(sampler, n = 100)
drawSamplesByStepCount <- function(sampler,
                                   n,
                                   burnin,
                                   burninleft,
                                   lengthout,
                                   thin,
                                   verbose){
  # Note: the initial graph is NOT returned at the moment
  samples <- vector("list", lengthout)
  if (verbose){
    cat("Drawing", n, "samples\n")
    flush.console()
    progress <- txtProgressBar(max = n, style = 3)
    setTxtProgressBar(progress, 0)
  }

  start <- proc.time()
  i <- 1
  while (i <= n){
    out <- sampler(i, burnin = burnin)

    if (i %% thin == 0 && i > burninleft){
      class(out) <- c("bn", "parental")
      samples[[(i-burninleft)/thin]] <- out
    }

    if (verbose){
      setTxtProgressBar(progress, i)
    }
    i <- i + 1
  }

  if (verbose){
    close(progress)
    elapsed <- (proc.time() - start)[[3]]
    cat("Drew", n, "samples in", elapsed, "seconds")
  }

  class(samples) <- c("mcmcbn", "bn.list", "parental.list")
  samples
}

#' List of MCMC Samplers
#'
#' method description
#'
#' @param ... Further arguments passed to method
#' @export
samplers <- function(...){
  x <- list(...)
  class(x) <- "samplers"
  x
}

#' Extract statistics from a sampler.
#'
#' Extracts the statistics collected during an MCMC run
#'
#' @param x A sampler
#' @param ... Further arguments passed to method
#' @return A list of statistics collected during the MCMC run.
#' @export
#' @seealso \code{\link{statistics.sampler}}
statistics <- function(x, ...){
  UseMethod("statistics")
}

#' Extract statistics from a sampler.
#'
#' Extracts the statistics collected during an MCMC run
#'
#' @param x A sampler
#' @param names Which statistics to extract? A character vector of
#'   statistics collected during the MCMC run.
#' @param ... Further arguments passed to method
#' @return A list of statistics collected during the MCMC run.
#' @export
#' @seealso \code{\link{statistics}}
statistics.sampler <- function(x, names, ...){
  stopifnot(inherits(x, "sampler"),
            length(names) > 0,
            inherits(names, "character"))

  statisticsTable <- get("statisticsTable", envir = environment(x))
  statisticsTable <- statisticsTable[, names, drop = F]
  nSteps <- get("nSteps", envir = environment(x))
  nBurnin <- get("nBurnin", envir = environment(x))
  statisticsTable[seq_len(nSteps - nBurnin), , drop = T]
}

#' Number of samples drawn.
#'
#' Returns the number of samples (MCMC steps) drawn in the supplied sampler.
#'
#' @param x A sampler
#' @param ... Further arguments, currently unused
#' @return The number of samples (steps) that have been drawn.
#' @export
#' @method length sampler
length.sampler <- function(x, ...){
  stopifnot(inherits(x, "sampler"))
  steps(x) - burnin(x)
}

#' Extract burnin
#'
#' Extracts the amount of burnin used in a sampler
#'
#' @param x A sampler
#' @param ... Further arguments passed to method
#' @return An integer amount of burnin
#' @export
#' @seealso \code{\link{steps.sampler}}
burnin <- function(x, ...){
  UseMethod("burnin")
}

#' Extract burnin
#'
#' Extracts the amount of burnin used in a sampler
#'
#' @param x A sampler
#' @param ... Further arguments, currently unused
#' @return An integer amount of burnin
#' @export
#' @method burnin sampler
burnin.sampler <- function(x, ...){
  stopifnot(inherits(x, "sampler"))
  get("nBurnin", envir = environment(x))
}

#' Extract burnin
#'
#' Extracts the amount of burnin used in a sampler
#'
#' @param x A sampler
#' @param ... Further arguments passed to method
#' @return An integer amount of burnin
#' @export
#' @seealso \code{\link{steps.sampler}}
burnin <- function(x, ...){
  UseMethod("burnin")
}

#' Number of samples drawn.
#'
#' Returns the number of samples (MCMC steps) drawn in the supplied sampler.
#'
#' @param x A sampler
#' @param ... Further arguments passed to method
#' @return The number of samples (steps) that have been drawn.
#' @export
#' @seealso \code{\link{steps.sampler}}
steps <- function(x, ...){
  UseMethod("steps")
}

#' Number of samples drawn.
#'
#' Returns the number of samples (MCMC steps) drawn in the supplied sampler.
#'
#' @param x A sampler
#' @param ... Further arguments, currently unused
#' @return The number of samples (steps) that have been drawn.
#' @export
#' @method length sampler
steps.sampler <- function(x, ...){
  stopifnot(inherits(x, "sampler"))
  get("nSteps", envir = environment(x))
}

#' Print a sampler
#'
#' Prints a sampler
#'
#' @param x A sampler
#' @param ... Further arguments, currently unused
#' @export
#' @method print sampler
print.sampler <- function(x, ...){
  cat("Number of steps: ", length(x))
  invisible(x)
}

#' Summarising MCMC samplers
#'
#' \code{summary} method for class "\code{sampler}"
#'
#' @param x A sampler
#' @param ... Further arguments, currently unused
#' @export
#' @method summary sampler
summary.sampler <- function(x, ...){
  cat("Number of steps: ", length(x))
  invisible(x)
}

#' Retreive MAP from sampler
#'
#' \code{map} method for class "\code{sampler}"
#'
#' @param x A sampler
#' @param ... Further arguments, currently unused
#' @export
#' @method map sampler
map.sampler <- function(x, ...){
  count <- get("count", envir = environment(x))
  count <- unlist(as.list(count))
  id <- names(count)[order(count, decreasing = T)[1]]
  as.bn(id)
}

#' Retreive top graph from sampler
#'
#' \code{top} method for class "\code{sampler}"
#'
#' @param x A sampler
#' @param head The number of graphs to consider
#' @param ... Further arguments, currently unused
#' @export
#' @method top sampler
top.sampler <- function(x, head = 10, ...){
  count <- get("count", envir = environment(x))
  count <- unlist(as.list(count))
  o <- order(count, decreasing = T)
  s <- seq_len(min(head, length(count)))
  id <- names(count)[o[s]]
  as.bn(id)
}

#' Retreive graph probabilities from sampler
#'
#' \code{gp} method for class "\code{sampler}"
#'
#' @param x A sampler
#' @param head The number of graphs to consider
#' @param ... Further arguments, currently unused
#' @export
#' @method gp sampler
gp.sampler <- function(x, head = 10000, ...){
  count <- get("count", envir = environment(x))
  count <- unlist(as.list(count))
  o <- order(count, decreasing = T)
  s <- seq_len(min(head, length(count)))
  count <- count[o[s]]
  count/sum(count)
}
