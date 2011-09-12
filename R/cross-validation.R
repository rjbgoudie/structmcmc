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
#' ..
#'
#' @param x A BN
#' @param train A data frame of training data
#' @param test A data frame of test data
#' @param metric A function that measures the distance between the
#'   predictions and the true values
#' @param verbose Logical indicating whether verbose output should be given
#' @export
residualsMultDir <- function(x, train, test, metric = kronecker_delta, verbose = F){
  stopifnot(inherits(x, "bn"),
            inherits(train, "data.frame"),
            inherits(test, "data.frame"))
  param <- bayes(x, train)
  sapply(seq_along(x), residualsMultDirNode,
                       x       = x,
                       param   = param,
                       test    = test,
                       metric  = metric,
                       verbose = verbose)
}

#' Residuals for a single node for a Multinomial-Dirichlet model
#' 
#' @inheritParams residualsMultDir
#' @param node An integer, giving the node
#' @export
residualsMultDirNode <- function(node, x, param, test, metric, verbose){
  pred <- predictNode(node, x, param, test)
  metric(pred, test[, node])
}

#' Predict a Multinomial-Dirichlet model
#'
#' ..
#'
#' @param x A BN
#' @param param Parameters for the model
#' @param newdata A data frame of new data
#' @export
predict.bn <- function(x, param, newdata){
  sapply(seq_along(x), predictNode,
                       x       = x,
                       param   = param,
                       newdata = newdata)
}

#' Predict a Multinomial-Dirichlet model
#' 
#' @inheritParams predict.bn
#' @param node An integer, giving the node
#' @export
predictNode <- function(node, x, param, newdata){
  param <- param[[node]]
  parents <- x[[node]]
  if (length(parents) > 0){
    ol <- levels(newdata[, node])
    newdata <- newdata[, parents, drop = F]
    if (require(plyr)){
      predictions <- adply(newdata, 1, makeRowPredictions,
                                       predictors = parents,
                                       param = param)
      factor(predictions[, ncol(predictions)], levels = ol)
    } else {
      stop("need plyr")
    }
  } else {
    ol <- levels(newdata[, node])
    n <- nrow(newdata)
    probabilities <- param[[1]]
    predictions <- ceiling(probabilities * n)

    # adjust for rounding errors
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

    # output final factor
    out <- c()
    for (i in seq_along(predictions)){
      out <- c(out, rep(names(predictions)[i], predictions[i]))
    }
    factor(out, levels = ol)
  }
}

#' Make row predictions
#' 
#' @param row A row
#' @param predictors
#' @param param 
makeRowPredictions <- function(row, predictors, param){
  id <- paste(unname(unlist(row)), collapse = ",")
  probabilities <- param[[id]]
  names(probabilities)[which.max(probabilities)]
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
