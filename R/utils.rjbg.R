# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Convert a data frame of factors to integers.
#'
#' Converts a \code{data.frame} that consists columns of factors into a
#' \code{data.frame} consisting of integers
#'
#' @param x an object of class \code{data.frame}
#' @param useLevelNames logical, indicating whether the labels of the
#'   levels should be converted to integers. This only makes sense if the
#'   levels are integers stored as characters. e.g. factor(c("3", "2", "3"))
#' @return The data.frame with columns converted to integers
#' @export
#' @seealso \code{\link{intAsFDF}}
fdfAsInt <- function(x, useLevelNames = T){
  stopifnot(
    class(x) == "data.frame",
    all(sapply(x, class) == "factor")
  )
  data.frame(lapply(x, function(col){
    if (useLevelNames){
      as.integer(levels(col))[as.integer(col)]
    }
    else {
      # the 1L needed to ensure that
      # the output is integer not double!
      as.integer(col) - 1L
    }
  }))
}

#' Convert a data frame of integers to factors.
#'
#' Converts a \code{data.frame} that consists columns of integers into a
#' \code{data.frame} consisting of factos
#'
#' @param x an object of class \code{data.frame}
#' @return The data.frame with columns converted to factors
#' @export
#' @seealso \code{\link{fdfAsInt}}
intAsFDF <- function(x){
  stopifnot(
    class(x) == "data.frame",
    all(sapply(x, class) %in% c("integer", "numeric"))
  )

  data.frame(lapply(x, as.factor))
}

#' Undocumented.
#'
#' method description
#'
#' @param a ...
#' @export
logsumexp <- function(x){
  # returns log(sum(exp(a)))
  # see sach email 19 Mar 2009
  .Call("logSumExp", as.numeric(x), PACKAGE = "structmcmc")
}

#' Multiple intersection function
#'
#' Performs intersection across all the supplied arguments.
#'
#' @param x A vector containing a sequence of items (conceptually with no
#'   duplicated values).
#' @param y A vector of the same \code{\link{mode}} as \code{x}.
#' @param ... Further vectors of the same mode as \code{x}.
#' @return A vector of the same \code{\link{mode}} as \code{x} and \code{y}.
#' @export
#' @examples
#' x <- c(1, 2, 4, 5)
#' y <- c(2, 3, 4)
#' intersection(x, y)
#' z <- c(5, 4, 3)
#' intersection(x, y, z)
intersection <- function(x, y, ...){
  if (missing(y)){
    .Internal(unique(x             = .Internal(unlist(x,
                                                      recursive = T,
                                                      use.names = F)),
                     incomparables = F,
                     fromLast      = F,
                     nmax          = NA))
  }
   else {
    if (missing(...)){
      # intersect2
      .Internal(unique(x             = y[.Call("fmatch", x, y, 0L, NULL,
                                         PACKAGE = "fastmatch")],
                       incomparables = F,
                       fromLast      = F,
                       nmax          = NA))
    } else {
      intersect2(x, intersection(y, ...))
    }
  }
}

#' Most-recently used stack.
#'
#' Stack that is aware of which items have been used most recently.
#'
#' @param size An integer. The stack size
#' @export
new_stack <- function(size = 1000L)
{
  keys_used <- vector("character", size)
  cache <- NULL
  cache_reset <- function(){
    cache <<- new.env(hash = TRUE, emptyenv(), size = size)
  }
  cache_get <- function(key){
    keys_used <<- keys_used[-match(key, keys_used)]
    keys_used <<- c(key, keys_used)
    get(key, env = cache, inherits = FALSE)
  }
  cache_has_key <- function(key){
    exists(key, env = cache, inherits = FALSE)
  }
  cache_rm <- function(key){
    if (cache_has_key(key)){
      rm(list = key, envir = cache)
    }
  }
  cache_set <- function(key, value){
    last <- keys_used[size]
    if (last != ""){
      cache_rm(last)
    }
    keys_used <<- keys_used[1:(size-1)]
    keys_used <<- c(key, keys_used)
    assign(key, value, env = cache)
  }
  cache_reset()
  list(reset = cache_reset,
       set = cache_set,
       get = cache_get,
       has_key = cache_has_key,
       keys = function() ls(cache),
       keys_used = function() keys_used)
}

#' Priority queue
#'
#' A priority queue, at least of sorts
#'
#' Based upon Hadley Wickham's memoise (MIT License)
#' http://cran.r-project.org/web/packages/memoise/
#'
#' @param threshold The amount below the top that should be kept
#' @export
new_queue <- function(threshold){
  queue <- list()
  priorities <- c()
  queue_reset <- function(){
    queue <<- vector("list", length = 0)
    priorities <<- vector("numeric", length = 0)
  }
  queue_list <- function(){
    queue
  }
  queue_which <- function(element){
    if (length(queue) > 0){
      which(sapply(queue, function(this){
        identical(this, element)
      }))
    } else {
      integer(0)
    }
  }
  queue_is_element <- function(element){
    length(queue_which(element)) > 0
  }
  queue_rm <- function(element){
    wh <- queue_which(element)
    if (length(wh) > 0){
      priorities <<- priorities[-wh]
      queue <<- queue[-wh]
    }
  }
  queue_top_priority <- function(element){
    max(priorities)
  }
  queue_prune <- function(){
    cutoff <- queue_top_priority() - threshold
    keep <- which(priorities > cutoff)
    queue <<- queue[keep]
    priorities <<- priorities[keep]
  }
  queue_add <- function(element, priority){
    if (!queue_is_element(element)){
      queue <<- c(queue, list(element))
      priorities <<- c(priorities, priority)
      queue_prune()
    }
  }
  queue_browser <- function(){
    browser()
  }
  queue_reset()
  list(reset = queue_reset,
       add = queue_add,
       ls = queue_list,
       rm = queue_rm,
       browser = queue_browser)
}
