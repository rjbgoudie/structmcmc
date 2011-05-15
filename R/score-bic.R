# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Undocumented.
#'
#' method description
#'
#' @param bn ...
#' @param data ...
#' @param nodes ...
#' @param mlcache ...
#' @param nijkcache ...
#' #' @export
scoreBIC <- function(bn,
                     data,
                     nodes = seq_along(bn),

                     mlcache = new.env(hash = T),
                     nijkcache = new.env(hash = T)){
  # add sanity checks

  # notation follows Heckerman (1999) in Learning in Graphical Models

  N <- nrow(data)

  getDim <- function(data, bn, nodes = seq_along(bn)){
    nl <- sapply(seq_len(ncol(data)), function(i) nlevels(data[, i]))

    getDimNode <- function(node){
      qi <- prod(nl[bn[[node]]])
      ri <- nl[node]
      qi * (ri - 1)
    }

    dims <- unlist(lapply(nodes, getDimNode))
    sum(dims)
  }
  d <- getDim(data, bn, nodes)

  theta <- unlist(ml(x              = bn,
                     data           = data,
                     nodes          = nodes,
                     regularisation = NaN,
                     cache          = mlcache))

  # fill in the blanks with rubbish.
  # the answer should be invariant to these, I think
  theta[is.nan(theta)] <- 0.5
  theta[theta == 0] <- 0.99999999

  Nijk <- unlist(getNijkCounts(data, bn, nodes, nijkcache))

  sum(Nijk * log(theta)) - d/2 * log(N)
}