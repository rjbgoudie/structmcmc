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
logsumexp <- function(a){
  # returns log(sum(exp(a)))
  # see sach email 19 Mar 2009

  m <- max(a)
  b <- a - m * rep(1, times = length(a))
  m + log(sum(exp(b)))
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
