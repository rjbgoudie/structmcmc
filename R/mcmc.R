# Part of the "structural" package, http://github.com/rjbgoudie/structural
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rjbgoudie/structural
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Draw samples from a MCMC sampler
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


#' Draw samples from a MCMC sampler, by time
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
    cat("Drawing samples, for ", time, " seconds")
    progress <- create_progress_bar("text")
    progress$init(time)
  }
  pb <- txtProgressBar(max = time, style = 3)

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
  close(pb)
  
  samples <- samples[seq_len(max(0,i-burninleft - 1))]
  class(samples) <- c("mcmcbn", "bn.list", "parental.list")
  samples
}

#' Draw samples from a MCMC sampler, by step count
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
    cat("Drawing", n, "samples")
    progress <- create_progress_bar("text")
    progress$init(n)
  }
  
  i <- 1
  while (i <= n){
    out <- sampler(i, burnin = burnin)

    if (i %% thin == 0 && i > burninleft){
      class(out) <- c("bn", "parental")
      samples[[(i-burninleft)/thin]] <- out
    }
    
    if (verbose){
      progress$step()
    }
    i <- i + 1
  }

  if (verbose){
    cat("\n")
  }

  class(samples) <- c("mcmcbn", "bn.list", "parental.list")
  samples
}
