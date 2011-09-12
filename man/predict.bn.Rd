\name{predict.bn}
\alias{predict.bn}
\title{Predict a Multinomial-Dirichlet model}
\usage{
  predict.bn(x, param, newdata, type = c("response",
  "probabilities"))
}
\arguments{
  \item{x}{A BN}

  \item{param}{Parameters for the model}

  \item{newdata}{A data frame of new data}

  \item{type}{The type of prediction required.}
}
\description{
  Predict each variable, given the model in the supplied
  Bayesian Network \code{x}.
}

