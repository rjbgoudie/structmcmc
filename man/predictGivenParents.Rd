\name{predictGivenParents}
\alias{predictGivenParents}
\title{Predict a node, given parents}
\usage{
  predictGivenParents(node, x, param, newdata, type)
}
\arguments{
  \item{node}{An integer, specifying the node}

  \item{x}{A BN}

  \item{param}{Parameters for the model}

  \item{newdata}{A data frame of new data}

  \item{type}{The type of prediction required.}
}
\description{
  Predicts a particular node, given the parents in a
  Bayesian network
}

