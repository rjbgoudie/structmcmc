# Part of the "structmcmc" package, http://github.com/rbtgde/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Multinomial-Dirichlet Log marginal likelihood.
#'
#' method description
#'
#' @export
#' @seealso \code{\link{logScoreZellnerFUN}}, \code{\link{logScoreMultDir}},
#'   \code{\link{logScoreZellner}}
logScoreMultDirFUN <- function(){
  list(offline = logScoreMultDirOffline,
       online  = logScoreMultDirIncremental,
       local   = localLogScoreMultDir,
       prepare = logScoreMultDirPrepare)
}

#' Normal Log marginal likelihood.
#'
#' method description
#'
#' @export
#' @seealso \code{\link{logScoreMultDirFUN}}, \code{\link{logScoreMultDir}},
#'   \code{\link{logScoreZellner}}
logScoreZellnerFUN <- function(){
  list(offline = logScoreZellnerOffline,
       online  = logScoreZellnerIncremental,
       local   = localLogScoreZellner,
       prepare = logScoreZellnerPrepare)
}