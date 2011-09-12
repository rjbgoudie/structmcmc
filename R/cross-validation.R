# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

#' Residuals from a Multinomial-Dirichlet model
#'
#' Given a Bayesian network, some training data and some test data,
#' the model given by fitting the Bayesian network to the training data is
#' used to predict each node of the test data, given the parents of that
#' node in the Bayesian network.
#' 
#' The residual is then computed, using the supplied \code{metric}.
#' 
#' Alternatively, a \code{bn.list} of Bayesian networks can be supplied,
#' together with a vector of weights. The models (Bayesian networks) are then
#' averaged over, according to the supplied weights, to give a model
#' averaging prediction.
#' 
#' The residuals are again computed, using the suppplied \code{metric}.
#'
#' @param x A BN, or a \code{bn.list}
#' @param weights A numeric vector weights for the models \code{x}
#' @param train A data frame of training data
#' @param test A data frame of test data
#' @param metric A function that measures the distance between the
#'   predictions and the true values
#' @param verbose Logical indicating whether verbose output should be given
#' @export
residualsMultDir <- function(x,
                             weights = 1,
                             train,
                             test,
                             metric = kronecker_delta,
                             verbose = F){
  is_bn <- inherits(x, "bn")
  stopifnot(is_bn || inherits(x, "bn.list"),
            inherits(train, "data.frame"),
            inherits(test, "data.frame"))
  if (is_bn){
    n <- length(x)
  } else {
    n <- length(x[[1]])
  }
  if (isTRUE(verbose)){
    progress <- txtProgressBar(max = n, style = 3)
    setTxtProgressBar(progress, 0)
  }
  for (node in seq_len(n)){
    out <- residualsMultDirNode(node    = node,
                                x       = x,
                                weights = weights,
                                train   = train,
                                test    = test,
                                metric  = metric,
                                verbose = verbose)
    if (isTRUE(verbose)){
      setTxtProgressBar(progress, node)
    }
  }
  if (isTRUE(verbose)){
    close(progress)
  }
  out
}

#' Residuals for a single node for a Multinomial-Dirichlet model
#'
#' Given a Bayesian network, some training data and some test data,
#' the model given by fitting the Bayesian network to the training data is
#' used to predict \code{node} of the test data, given the parents of that
#' node in the Bayesian network.
#' 
#' The residual is then computed, using the supplied \code{metric}.
#' 
#' Alternatively, a \code{bn.list} of Bayesian networks can be supplied,
#' together with a vector of weights. The models (Bayesian networks) are then
#' averaged over, according to the supplied weights, to give a model
#' averaging prediction. Only a single node is considered.
#' 
#' The residuals are again computed, using the suppplied \code{metric}.
#'
#' @inheritParams residualsMultDir
#' @param node An integer, giving the node
#' @export
residualsMultDirNode <- function(node,
                                 weights = 1,
                                 x,
                                 train,
                                 test,
                                 metric,
                                 verbose){
  if (inherits(x, "bn")){
    param <- bayes(x, train)
    pred <- predictNode(node = node,
                        x = x,
                        param = param,
                        test = test)
    metric(pred, test[, node])
  } else {
    pred <- predictModelAverageNode(node = node,
                                    x = x,
                                    weights = weights,
                                    train = train,
                                    newdata = test)
    metric(pred, test[, node])
  }
}

#' Predict a Multinomial-Dirichlet model
#'
#' Predict each variable, given the model in the supplied Bayesian Network
#' \code{x}.
#'
#' @param x A BN
#' @param param Parameters for the model
#' @param newdata A data frame of new data
#' @param type The type of prediction required.
#' @export
predict.bn <- function(x,
                       param,
                       newdata,
                       type = c("response", "probabilities")){
  sapply(seq_along(x), predictNode,
                       x       = x,
                       param   = param,
                       newdata = newdata)
}

#' Predict a Multinomial-Dirichlet model
#'
#' Predict a node of a Bayesian network, given a collection of Bayesian
#' networks. The predictions from each supplied Bayesian network are averaged
#' over, according to their weights.
#'
#' @inheritParams predict.bn
#' @param node An integer, giving the node
#' @export
predictModelAverageNode <- function(node,
                                    x,
                                    weights,
                                    train,
                                    newdata){
  probabilities <- lapply(x, function(model){
    param <- bayes(model, train)
    predictNode(node, model, param, newdata, type = "probabilities")
  })
  probabilities <- Map("*", probabilities, weights)
  probabilities <- Reduce("+", probabilities)
  response <- adply(probabilities, 1, rowProbabiltiesToResponse,
                                      names = levels(newdata[, node]))
  ol <- levels(newdata[, node])
  factor(response[, ncol(response)], levels = ol)
}

#' Get prediction for a row
#'
#' @param prob A numeric vector of probabilities
#' @param names A vector of level names
rowProbabiltiesToResponse <- function(prob, names){
  names[which.max(prob)]
}

#' Predict a Multinomial-Dirichlet model
#'
#' @inheritParams predict.bn
#' @param node An integer, giving the node
#' @export
predictNode <- function(node,
                        x,
                        param,
                        newdata,
                        type = c("response", "probabilities")){
  if (length(x[[node]]) > 0){
    predictGivenParents(node, x, param, newdata, type = type)
  } else {
    predictGivenNoParents(node, x, param, newdata, type = type)
  }
}

#' Adjust the predictions
#'
#' Due to rounding, too many or too few predictions may be made. This
#' functions makes the adjustments in a reasonable manner
#'
#' @param predictions An integer vector of the number of each category that
#'   have been predicted.
#' @param probabilities A vector of probabilities
#' @param n How many predictions should there be
adjustPredictions <- function(predictions, probabilities, n){
  while (sum(predictions) > n){
    o <- sort(decim(probabilities), decreasing = F)

    done <- F
    allowed <- seq_along(probabilities)
    while (!done){
      if (predictions[allowed[1]] > 0){
        predictions[allowed[1]] <- predictions[allowed[1]] - 1
        done <- T
      } else {
        allowed <- allowed[-1]
      }
    }
  }
  while (sum(predictions) < n){
    w <- which.max(decim(probabilities))
    predictions[w] <- predictions[w] + 1
  }
  predictions
}

#' Predict a node, given parents
#'
#' Predicts a particular node, given the parents in a Bayesian network
#'
#' @inheritParams predict.bn
#' @param node An integer, specifying the node
predictGivenParents <- function(node, x, param, newdata, type){
  param <- param[[node]]
  parents <- x[[node]]
  if (require(plyr)){
    predictions <- adply(newdata, 1, makeRowPredictions,
                                     predictors = parents,
                                     param = param,
                                     type = type)
    if (type[1] == "response"){
      ol <- levels(newdata[, node])
      factor(predictions[, ncol(predictions)], levels = ol)
    } else {
      predictions[-seq_len(ncol(newdata))]
    }
  } else {
    stop("need plyr")
  }
}

#' Predict a node, given no parents
#'
#' Predicts a particular node
#'
#' @inheritParams predict.bn
#' @param node An integer, specifying the node
predictGivenNoParents <- function(node,
                             x,
                             param,
                             newdata,
                             type = c("response", "probabilities")){
  param <- param[[node]]
  ol <- levels(newdata[, node])
  n <- nrow(newdata)
  probabilities <- param[[1]]

  if (type[1] == "response"){
    predictions <- ceiling(probabilities * n)

    # adjust for rounding errors
    predictions <- adjustPredictions(predictions, probabilities, n)

    out <- c()
    for (i in seq_along(predictions)){
      out <- c(out, rep(names(predictions)[i], predictions[i]))
    }
    factor(out, levels = ol)
  } else {
    matrix(rep(probabilities, each = nrow(newdata)), nrow = nrow(newdata))
  }
}

#' Make row predictions
#'
#' @param row A vector containing the predictors in the model
#' @param predictors The parents of the node
#' @param param A list of parameters
#' @param type The type of prediction required.
makeRowPredictions <- function(row,
                               predictors,
                               param,
                               type = c("response", "probabilities")){
  id <- paste(unname(unlist(row[predictors])), collapse = ",")
  probabilities <- param[[id]]
  if (type[1] == "response"){
    names(probabilities)[which.max(probabilities)]
  } else {
    probabilities
  }
}

#' Find the non-integer part of a number
#'
#' @param x A numeric vector
#' @examples
#' decim(2.3)
#' decim(2)
decim <- function(x){
  x - floor(x)
}

#' Residuals from a Normal model
#'
#' ...
#'
#' @param x A BN
#' @param train A data frame of training data
#' @param test A data frame of test data
#' @param metric A function that measures the distance between the
#'   predictions and the true values
#' @export
residsNormal <- function(x, train, test, metric = square_error){
  sapply(seq_along(x), residsNormalNode,
                       x, train, test, metric, verbose)
}

#' Residuals for a single node for a Normal model
#'
#' @inheritParams residsNormal
#' @param node An integer, giving the node
#' @export
residsNormalNode <- function(node, x, train, test, metric, verbose){
  parents <- x[[node]]
  if (length(parents) > 0){
    outcome <- colnames(dat)[node]
    predictors <- colnames(dat)[parents]

    outcomeQuote <- paste("`", outcome, "`", sep = "")
    predictorsQuote <- paste("`", predictors, "`", sep = "")
    predictorsStr <- paste(predictorsQuote, collapse = " + ")
    formula <- paste(outcomeQuote, " ~ ", predictorsStr)
    formula <- as.formula(formula)

    thisTraining <- train[, c(outcome, predictors), drop = F]
    fit <- zlm(formula, data = thisTraining)

    thisTest <- test[, predictors, drop = F]
    pred <- predict(fit, newdata = thisTest)

    metric(pred, test[, node])
  } else {
    NULL
  }
}

#' Metric: Kronecker delta
#'
#' @param predictions A factor variable
#' @param true A factor variable
#' @export
kronecker_delta <- function(predictions, true){
  stopifnot(inherits(predictions, "factor"),
            inherits(true, "factor"),
            identical(levels(predictions), levels(true)),
            identical(length(predictions), length(true)))
  sum(predictions != true)
}

#' Metric: Square error
#' 
#' @param predictions A numeric variable
#' @param true A numeric variable
#' @export
square_error <- function(predictions, true){
  differences <- predictions - true
  sum(differences^2)
}
