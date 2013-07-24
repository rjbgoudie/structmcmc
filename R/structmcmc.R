# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
#
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
#
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' The structmcmc package.
#'
#' @import parental lattice fastmatch bit
#' @docType package
#' @name structmcmc
#' @aliases structmcmc package-structmcmc
#' @useDynLib structmcmc
#' @seealso \code{\link{posterior}}
#' @examples
#' # Setup data frame
#' x1 <- factor(c("a", "a", "g", "c", "c", "a", "g", "a", "a"))
#' x2 <- factor(c(2, 2, 4, 3, 1, 4, 4, 4, 1))
#' x3 <- factor(c(FALSE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
#' x <- data.frame(x1 = x1, x2 = x2, x3 = x3)
#'
#' # Draw samples from the posterior using MC^3.
#' set.seed(1234)
#' mcmc <- posterior(data = x, method = "mh-mcmc",
#'                   nSamples = 1000, nBurnin = 100)
#'
#' # Compute and plot estimated edge probabilities
#' epmcmc <- ep(mcmc)
#'
#' # Exact evaluation by exhaustive enumeration
#' exact <- posterior(x, "exact")
#' epexact <- ep(exact)
#'
#' # Comparing multiple MCMC runs
#' mcmc2 <- posterior(x, "mh-mcmc", nSamples = 1000, nBurnin = 100)
#' epmcmc2 <- ep(mcmc2)
NULL

#' Internal function.
#'
#' @param libname ...
#' @param pkgname ...
#' @name internal-onload
#' @aliases internal-onload .onLoad
.onLoad <- function(libname, pkgname){
  cat("copyright (c) 2008, Robert J. B. Goudie, University of Warwick\n")
  cat('For citation information, type citation("structmcmc").\n')
  cat('Type help("package-structmcmc") to get started.\n')
}
