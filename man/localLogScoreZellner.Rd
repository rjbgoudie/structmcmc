\name{localLogScoreZellner}
\alias{localLogScoreZellner}
\title{Local Normal-inverse-gamma (with g-prior) Log marginal likelihood.}
\usage{localLogScoreZellner(node, parents, logScoreParameters,
    cache, checkInput=T)}

\description{
  Local Normal-inverse-gamma (with g-prior) Log marginal
  likelihood.
}

\details{
  Compute the LOCAL log marginal likelihood of the supplied
  Bayesian Networks. ie the contribution to the log
  marginal liklihood from one individual node.

  Let \eqn{X}{X} be a data matrix with a number of
  predictors (in columns), and \eqn{y}{y} be an response
  variable, and that \eqn{n}{n} observations are available
  for each. For a graph \eqn{G}{G} (since this is local
  score this is equivalent to an indicator vector), the
  model used is takes the form \eqn{y = \phi_{G} \beta +
  \epsilon}{y = phi_G * beta + epsilon} with \eqn{\epsilon
  \sim N(0, \sigma^{2} I)}{epsilon ~ N(0, sigma^{2} I)}.
  Note that the data needs to be standardised
  (zero-meaned).

  The design matrix \eqn{\phi_{G}}{phi_{G}} is a column of
  1s, and then columns corresponding to each of the parents
  of the node. No cross-terms are included.

  The prior used factorises as \eqn{p(\beta, \sigma) =
  p(\beta \mid \sigma)p(\sigma)}{ p(beta, sigma) = p(beta |
  sigma)p(sigma)}, The variance has an uninformative, scale
  invariant Jeffrey's prior \eqn{p(\sigma) =
  1/\sigma^{2}}{p(sigma) = 1/sigma^2}, and the coefficients
  have a zero-mean Normal prior (a Zellner g-prior), with
  \eqn{g = n}, so that \eqn{\beta \mid \sigma \sim N(0, g
  \sigma^{2}\left(\phi_{G}^{\prime}\phi_{G}\right)^{-1})}{
  beta | sigma ~ N(0, g * sigma^2 * (phi'_G phi_G)^-1)}

  The above specification gives the following marginal
  likelihood.

  \deqn{p(y \mid G) \propto (1 + n)^{-(\eta + 1)/2}
  \left(X^{T} X - \frac{n}{n + 1} X^{T}
  \phi_{G}(\phi^{T}_{G}\phi_{G})^{-1}\phi^{T}_{G}
  X\right)^{-\frac{n}{2}} }{ P(y | G) propto (1 + n)^(-(eta
  + 1)/2) * (X' * X - (n/(n + 1)) * X' * phi_G * (phi'_G *
  phi_G)^(-1) * phi_G * X)^(-n/2) }
}
\value{A numeric vector of length 1, giving the log marginal likelihood.
The environment 'cache' will also be updated because its scope is
global.}
\references{Nott, D. J., & Green, P. J. (2004). \emph{Bayesian Variable Selection
and the Swendsen-Wang Algorithm}. Journal of Computational and Graphical
Statistics, 13, 141-157.
\url{http://dx.doi.org/10.1198/1061860042958}

Smith, M., & Kohn, R. (1996). \emph{Nonparametric Regression using
Bayesian Variable Selection}. Journal of Econometrics, 75, 317-343.
\url{http://dx.doi.org/10.1016/0304-4076(95)01763-1}.

Geiger, D., & Heckerman, D. (1994). \emph{Learning Gaussian Networks}.
Proceedings of the 10th Conference Annual Conference on Uncertainty in
Artificial Intelligence (UAI-94), 235-240.
\url{http://uai.sis.pitt.edu/displayArticleDetails.jsp?mmnu=1&smnu=2&
article_id=509&proceeding_id=10}}
\seealso{\code{\link{logScoreZellner}},
\code{\link{logScoreZellnerOffline}},
\code{\link{logScoreZellnerIncremental}}.}
\arguments{
  \item{node}{A numeric vector of length 1. The node to compute the local
log score for.}
  \item{parents}{A numeric vector. The parents of node.}
  \item{logScoreParameters}{A list with the following components:
\describe{
\item{data}{A matrix with columns giving the values of each random
variable.}
\item{nl}{A numeric vector of length nNodes(currentBN), specifying the
number of levels that each random variable takes.}
}}
  \item{cache}{Optionally, provide an environment with cached local scores
for this data.}
  \item{checkInput}{A logical of length 1, specifying whether to check the
inputs to the function.}
}
