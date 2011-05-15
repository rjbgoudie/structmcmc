# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Convert to grain.
#'
#' This is an interface to the package \pkg{gRain}. We supply a Bayesian
#' Network, its parameters, and the raw data. An object of class 'grain'
#' is returned. This object can be queried
#'
#' Specifically, this function is a wrapper around
#' \code{\link[gRain]{cptable}}, \code{\link[gRain]{extractCPT}} and
#' \code{\link[gRain]{grain}}.
#'
#' @param x A list of multinomial parameters, in the form
#'   output by \code{\link{ml}}
#' @param net a Bayesian Network. An object of class 'bn'
#' @param dat A data.frame with factors in columns. The columns should
#'   correspond to the nodes of net, and the parameters of x.
#' @return An object of class "grain".
#' @seealso \code{\link[gRain]{querygrain}} for querying the result. For
#'   queries involving conditioning, use \code{\link[gRain]{setFinding}}
#'   before running \code{querygrain}.
#' @export
#' @seealso \code{\link{marginalGivenOthers}},
#'   \code{\link{marginalGivenIntervention}}, \code{\link{queryFinding}}
as.grain <- function(x, net, dat){
  stopifnot(class(x)   ==   "list",
            "bn"       %in% class(net),
            class(dat) ==   "data.frame")
  if (require(gRain)){
    nNodes <- nNodes(net)
    nodeNames <- names(dat)
    ll <- lapply(dat, levels)
    cptableList <- vector("list", length = nNodes)
    for (i in 1:nNodes){
      thisNodeName <- nodeNames[i]
      parents <- net[[i]]
      nParents <- length(parents)
      parentNames <- nodeNames[parents]
      # construct the formula for cptable
      # this specifies the node, and its parents
      parentsString <- paste("`", parentNames, "`",
                             collapse = " + ", sep = "")
      expr <- paste("`", thisNodeName, "`", sep = "")
      if (nParents > 0){
        expr <- paste(expr, parentsString, sep = " | ")
      }
      expr <- paste("~", expr)
      form <- as.formula(expr)
      values <- unlist(x[[i]])
      # NaN is output by ml() when there are no observations
      # for that parent set.
      if (any(is.nan(values))){
        warning("Supplied parameter set 'x' contains NaN")
        values[is.nan(values)] <- 1
      }
      cptableList[[i]] <- cptable(form,
                                  values = values,
                                  levels = ll[[i]])
    }
    plist <- compileCPT(cptableList)
    grain(plist)
  } else {
    stop("The package gRain must be installed")
  }
}

#' Get all marginal probabilities.
#'
#' Get the marginal probability of a node, when each other node (separately)
#' has been conditioned upon each of its levels.
#'
#' @param x A grain object
#' @param node An integer or a name of a column of \code{dat}
#' @param dat A data frame, with columns corresponding to the Bayes Net
#'   in \code{grain}.
#' @param FUN a function that is applied to the output of
#'   \code{\link[gRain]{querygrain}}. The first argument should accept
#'   the result of \code{\link[gRain]{querygrain}} FOR A PARTICULAR NODE.
#'   (Usaully \code{querygrain} returns a list, since it may provide the
#'   distribution of a number of nodes at once. Here, we remove that list.)
#'   Any other arguments are passed \code{...}
#' @param ... Passed to \code{FUN}
#' @return A list, with a component for each node of the Bayes Net (apart
#'   from the node \code{node}). Each component is a list, corresponding to
#'   the named level. These components contain the result of applying
#'   result of applying \code{FUN} to the result of \code{querygrain}.
#' @export
#' @seealso \code{\link{marginalGivenIntervention}},
#'   \code{\link{as.grain}}, \code{\link{queryFinding}}
marginalGivenOthers <- function(x, node, dat, FUN = identity, ...){
  stopifnot("grain"      %in% class(x),
            class(node)  %in% c("numeric", "character"),
            length(node) ==   1,
            class(dat)   ==   "data.frame",
            class(FUN)   ==   "function")
  if (class(node) == "numeric"){
    node <- names(dat)[node]
  }
  nodeID <- which(names(dat) == node)
  ll <- lapply(dat, levels)[-nodeID]
  nn <- names(dat)[-nodeID]
  names(ll) <- nn
  out <- lapply(nn, function(nodeWithState){
    thisNodeLevels <- ll[[nodeWithState]]
    result <- lapply(thisNodeLevels, function(nodeState){
      queryFinding(x             = x,
                   node          = node,
                   nodeWithState = nodeWithState,
                   nodeState     = nodeState,
                   FUN           = FUN)
    })
    names(result) <- thisNodeLevels
    result
  })
  names(out) <- nn
  out
}

#' Query an independence network, given a finding.
#'
#' Query a finding, and apply \code{FUN}. Just a wrapper around
#' \code{\link[gRain]{setFinding}} and \code{\link[gRain]{querygrain}}.
#'
#' @param x A grain
#' @param node Character vector of length 1. The node whose marginal
#'   distribution is required.
#' @param nodeWithState Character vector of length 1. The node whose state
#'   is to be conditioned upon.
#' @param nodeState The state to set \code{nodeWithState} to. A character
#'   vetor of length 1.
#' @param FUN a function that is applied to the output of
#'   \code{\link[gRain]{querygrain}}. The first argument should accept
#'   the result of \code{\link[gRain]{querygrain}} FOR A PARTICULAR NODE.
#'   (Usaully \code{querygrain} returns a list, since it may provide the
#'   distribution of a number of nodes at once. Here, we remove that list.)
#'   Any other arguments are passed \code{...}
#' @param ... Passed to \code{FUN}
#' @return The result of \code{FUN}
#' @export
#' @seealso \code{\link{as.grain}}, \code{\link{marginalGivenOthers}},
#'   \code{\link{marginalGivenIntervention}}
queryFinding <- function(x, node, nodeWithState, nodeState, FUN, ...){
  stopifnot("grain"               %in% class(x),
            class(node)           ==   "character",
            length(node)          ==   1,
            class(nodeWithState)  ==   "character",
            length(nodeWithState) ==   1,
            class(nodeState)      ==   "character",
            class(FUN)            ==   "function")
  x <- setFinding(object = x,
                  nodes  = nodeWithState,
                  states = nodeState)
  tab <- querygrain(object = x,
                    nodes  = node,
                    type   = "marginal")[[1]]
  FUN(tab, ...)
}

#' Get marginal probabilities, given an intervention.
#'
#' Get the marginal probability of a node, when each other node (separately)
#' has been conditioned upon each of its levels.
#'
#' @param net A BN. A "parental" object.
#' @param node An integer or a name of a column of \code{dat}
#' @param dat A data frame, with columns corresponding to the Bayes Net
#'   in \code{grain}.
#' @param FUN a function that is applied to the output of
#'   \code{\link[gRain]{querygrain}}. The first argument should accept
#'   the result of \code{\link[gRain]{querygrain}} FOR A PARTICULAR NODE.
#'   (Usaully \code{querygrain} returns a list, since it may provide the
#'   distribution of a number of nodes at once. Here, we remove that list.)
#'   Any other arguments are passed \code{...}
#' @param ... Passed to \code{FUN}
#' @return A list, with a component for each node of the Bayes Net (apart
#'   from the node \code{node}). Each component is a list, corresponding to
#'   the named level. These components contain the result of applying
#'   result of applying \code{FUN} to the result of \code{querygrain}.
#' @export
#' @seealso \code{\link{marginalGivenOthers}}, \code{\link{as.grain}},
#'   \code{\link{queryFinding}}
marginalGivenIntervention <- function(net,
                                      node,
                                      dat,
                                      FUN = identity,
                                      ...){
  stopifnot("bn"         %in% class(net),
            class(node)  %in% c("numeric", "character"),
            length(node) ==   1,
            class(dat)   ==   "data.frame",
            class(FUN)   ==   "function")
  if (class(node) == "numeric"){
    node <- names(dat)[node]
  }
  nodeID <- which(names(dat) == node)
  ll <- lapply(dat, levels)[-nodeID]
  nn <- names(dat)[-nodeID]
  names(ll) <- nn
  out <- lapply(seq_along(nn), function(i){
    nodeWithState <- nn[i]
    net[[i]] <- integer(0)
    estimate <- bayes(net, dat, prior = "qi")
    x <- as.grain(estimate, net, dat)
    
    thisNodeLevels <- ll[[nodeWithState]]
    result <- lapply(thisNodeLevels, function(nodeState){
      queryFinding(x             = x,
                   node          = node,
                   nodeWithState = nodeWithState,
                   nodeState     = nodeState,
                   FUN           = FUN)
    })
    names(result) <- thisNodeLevels
    result
  })
  names(out) <- nn
  out
}
