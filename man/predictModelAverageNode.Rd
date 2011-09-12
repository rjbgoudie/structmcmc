\name{predictModelAverageNode}
\alias{predictModelAverageNode}
\title{Predict a Multinomial-Dirichlet model}
\usage{
  predictModelAverageNode(node, x, weights, train, newdata)
}
\arguments{
  \item{node}{An integer, giving the node}

  \item{x}{A BN}

  \item{newdata}{A data frame of new data}
}
\description{
  Predict a node of a Bayesian network, given a collection
  of Bayesian networks. The predictions from each supplied
  Bayesian network are averaged over, according to their
  weights.
}

