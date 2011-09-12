\name{residualsMultDirNode}
\alias{residualsMultDirNode}
\title{Residuals for a single node for a Multinomial-Dirichlet model}
\usage{
  residualsMultDirNode(node, weights = 1, x, train, test,
  metric, verbose)
}
\arguments{
  \item{node}{An integer, giving the node}

  \item{weights}{A numeric vector weights for the models
  \code{x}}

  \item{x}{A BN, or a \code{bn.list}}

  \item{train}{A data frame of training data}

  \item{test}{A data frame of test data}

  \item{metric}{A function that measures the distance
  between the predictions and the true values}

  \item{verbose}{Logical indicating whether verbose output
  should be given}
}
\description{
  Given a Bayesian network, some training data and some
  test data, the model given by fitting the Bayesian
  network to the training data is used to predict
  \code{node} of the test data, given the parents of that
  node in the Bayesian network.
}
\details{
  The residual is then computed, using the supplied
  \code{metric}.

  Alternatively, a \code{bn.list} of Bayesian networks can
  be supplied, together with a vector of weights. The
  models (Bayesian networks) are then averaged over,
  according to the supplied weights, to give a model
  averaging prediction. Only a single node is considered.

  The residuals are again computed, using the suppplied
  \code{metric}.
}

