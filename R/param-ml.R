# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Maximum likelihood estimates.
#'
#' A generic
#'
#' @param x An object
#' @param ... Further arguments passed to method
#' @export
ml <- function (x, ...) {
  UseMethod("ml")
}

#' Maximum likelihood estimates for parameters of a BN.
#'
#' Computes the maximum likelihood estimates for the parameters of a
#' Bayesian Network. These are just the proportions of each category for
#' each configuration of the parents of a node in the Bayesian Network.
#'
#' If, for a particular configuration of the parents, we have no
#' observations, then by default, \code{NaN} is returned.
#'
#' @param x The Bayesian Network. An object of class 'bn'
#' @param data A data frame
#' @param nodes A subset of 1, ..., \code{nNodes(x)}. A numeric vector.
#' @param regularisation One of \code{NaN}, \code{"qi"}, or a numeric vector
#'   of length 1. Supplying \code{NaN} will place a \code{NaN} in any
#'   parts of the table where there is no information. Supplying
#'   \code{"qi"} will add a factor of \code{1/(number of parents of node)} to
#'   each value (as in Bayesian inference). Supplying a number will add this
#'   to each value.
#' @param cache A cache
#' @param ... Further arguments (unused)
#' @return A list of length \code{nNodes(x)}. Each component is a list
#'   containing components for each configuration of that node's parents
#'   in the Bayesian Network \code{x}. Each of these components is a
#'   numeric vector of probabilities that sum to 1, labelled with the
#'   levels of the relevant node.
#' @S3method ml bn
#' @method ml bn
ml.bn <- function(x,
                  data,
                  nodes = seq_along(x),
                  regularisation = NaN,
                  cache = new.env(hash = T),
                  ...){
  getNodeMLParameter <- function(data, node, parents, regularisation, cache){
    id <- .Internal(paste(list(c(node, parents)), "", ","))
    cacheRecord <- cache[[id]]
    if (!is.null(cacheRecord)){
      cacheRecord
    }
    else {
      nl <- sapply(data, nlevels)

      dataSubset <- data[, c(node, parents)]
      dataSubsetNames <- names(dataSubset)
      rawCounts_ijk <- table(dataSubset, dnn = dataSubsetNames)

      if (!is.nan(regularisation)){
        if (regularisation == "qi"){
          prior <- 1/prod(nl[parents])
          rawCounts_ijk <- rawCounts_ijk + prior
        } else {
          rawCounts_ijk <- rawCounts_ijk + regularisation
        }
      }

      parentsSeq <- seq_along(parents) + 1
      nodeNLevels <- nlevels(data[, node])

      pt <- prop.table(rawCounts_ijk, margin = parentsSeq)
      parentsNLevels <- prod(dim(pt))/nodeNLevels

      dim(pt) <- c(nodeNLevels, parentsNLevels)
      rownames(pt) <- levels(data[, node])
      result <- unlist(apply(pt, 2, list), rec = F)

      # deal with names
      # if no parents, don't do any naming
      if (length(parents) == 1){
        names(result) <- levels(data[, parents])
      }
      else if (length(parents) > 1){
        parentLevelNamesList <- lapply(parents, function(i) levels(data[,i]))

        egrid <- do.call("expand.grid", parentLevelNamesList)
        dimNames <- apply(egrid, 1, paste, collapse = ",")

        names(result) <- dimNames
      }

      (cache[[id]] <- result)
    }
  }

  result <- lapply(nodes, function(node){
    getNodeMLParameter(data, node, x[[node]], regularisation, cache)
  })
  names(result) <- nodes
  result
}

#' Undocumented.
#'
#' method description
#'
#' @param data ...
#' @param bn ...
#' @param nodes ...
#' @param cache ...
#' @export
getNijkCounts <- function(data, bn, nodes = seq_along(bn),
                          cache = new.env(hash = T)){
  getNodeNijkCounts <- function(data, node, parents, cache){
    id <- .Internal(paste(list(c(node, parents)), "", ","))
    cacheRecord <- cache[[id]]
    if (!is.null(cacheRecord)){
      cacheRecord
    }
    else {
      dataSubset <- data[, c(node, parents)]
      dataSubsetNames <- names(dataSubset)
      rawCounts_ijk <- table(dataSubset, dnn = dataSubsetNames)

      parentsSeq <- seq_along(parents) + 1
      nodeNLevels <- nlevels(data[, node])

      parentsNLevels <- prod(dim(rawCounts_ijk))/nodeNLevels
      dim(rawCounts_ijk) <- c(nodeNLevels, parentsNLevels)
      rownames(rawCounts_ijk) <- levels(data[, node])

      result <- unlist(apply(rawCounts_ijk, 2, list), rec = F)

      # deal with names
      # if no parents, don't do any naming
      if (length(parents) == 1){
        names(result) <- levels(data[, parents])
      }
      else if (length(parents) > 1) {
        parentLevelNamesList <- lapply(dataSubset[, -1], function(i){
          list(levels(i))
        })

        # the list's names confuse do.call
        # it thinks they are argument names
        parentLevelNamesList <- unname(parentLevelNamesList)

        egrid <- do.call("expand.grid", parentLevelNamesList)
        dimNames <- apply(egrid, 1, paste, collapse = ",")

        names(result) <- dimNames
      }
      (cache[[id]] <- result)
    }
  }

  result <- lapply(nodes, function(node){
    getNodeNijkCounts(data, node, bn[[node]], cache)
  })
  names(result) <- nodes
  result
}
