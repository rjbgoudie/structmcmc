\name{predictNode}
\alias{predictNode}
\title{Predict a Multinomial-Dirichlet model}
\usage{
  predictNode(node, x, param, newdata,
    type = c("response", "probabilities"))
}
\arguments{
  \item{node}{An integer, giving the node}

  \item{x}{A BN}

  \item{param}{Parameters for the model}

  \item{newdata}{A data frame of new data}

  \item{type}{The type of prediction required.}
}
\description{
  Predict a Multinomial-Dirichlet model
}

