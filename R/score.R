# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Multinomial-Dirichlet Log marginal likelihood.
#'
#' This function returns a list of functions that are used for MCMC
#' computation for Multinomial-Dirichlet models.
#'
#' @return A list consisting of the functions that perform each of the
#' following roles.
#' \describe{
#'   \item{offline}{A function that computes the logScore of a Bayesian
#'                  Network}
#'   \item{online}{A function that incrementally computes the logScore of a
#'                 Bayesian Network}
#'   \item{local}{A function that computes the local logScore of a
#'                Bayesian Network}
#'   \item{prepare}{A function that prepares the data, and any further
#'                  pre-computation required by the logScore functions.}
#' }
#'
#' @export
#' @seealso \code{\link{logScoreNormalFUN}}, \code{\link{logScoreMultDir}},
#'   \code{\link{logScoreNormal}}
logScoreMultDirFUN <- function(){
  list(offline = logScoreMultDirOffline,
       online  = logScoreMultDirIncremental,
       local   = localLogScoreMultDir,
       prepare = logScoreMultDirPrepare)
}

#' Normal Log marginal likelihood.
#'
#' This function returns a list of functions that are used for MCMC
#' computation for Normal models, with Zellner g-priors.
#'
#' @return A list consisting of the functions that perform each of the
#' following roles.
#' \describe{
#'   \item{offline}{A function that computes the logScore of a Bayesian
#'                  Network}
#'   \item{online}{A function that incrementally computes the logScore of a
#'                 Bayesian Network}
#'   \item{local}{A function that computes the local logScore of a
#'                Bayesian Network}
#'   \item{prepare}{A function that prepares the data, and any further
#'                  pre-computation required by the logScore functions.}
#' }
#'
#' @export
#' @seealso \code{\link{logScoreMultDirFUN}}, \code{\link{logScoreMultDir}},
#'   \code{\link{logScoreNormal}}
logScoreNormalFUN <- function(){
  list(offline = logScoreNormalOffline,
       online  = logScoreNormalIncremental,
       local   = localLogScoreNormal,
       prepare = logScoreNormalPrepare)
}
