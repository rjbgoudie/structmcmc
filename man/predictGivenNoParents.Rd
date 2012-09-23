\name{predictGivenNoParents}
\alias{predictGivenNoParents}
\title{Predict a node, given no parents}
\usage{
  predictGivenNoParents(node, x, param, newdata,
    type = c("response", "probabilities"))
}
\arguments{
  \item{node}{An integer, specifying the node}

  \item{x}{A BN}

  \item{param}{Parameters for the model}

  \item{newdata}{A data frame of new data}

  \item{type}{The type of prediction required.}
}
\description{
  Predicts a particular node
}

