# Part of the "structmcmc" package, https://github.com/rjbgoudie/structmcmc
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   https://github.com/rjbgoudie/structmcmc
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

# only handles binary???
# maybe only binary in the response?

#' Undocumented.
#'
#' method description
#'
#' @param graph ...
#' @param response ...
#' @param data ...
#' @param fold ...
#' @param verbose ...
#' @param binomial ...
#' @param predictiveReplications ...
#' @param ... Further arguments passed to method
#' @param Freq This argument was missing, and was picked up by R-check
#' @export
dagcv <- function(graph, response, data, fold = 10, verbose = F,
                  binomial = F, predictiveReplications = 10000, Freq, ...){
  stopifnot(
    class(data) == "data.frame",
    all(unlist(lapply(data, class)) == "factor"),
    "bn" %in% class(graph)
  )
  
  if (require(nnet) && require(reshape)){
    numberOfCases <- nrow(data)

    if (is.list(fold)){
      foldSeq <- seq_along(fold)
      folds <- fold
      fold <- length(fold)
    }
    else if (is.numeric(fold) && length(fold) == 1){
      foldSeq <- seq_len(fold)
      ordering <- sample.int(numberOfCases, numberOfCases, replace = F)
      folds <- split(ordering, foldSeq)
    }
    else {
      stop("fold not of recognised type")
    }
    foldLengths <- unlist(lapply(folds, length))

    numberOfTrials <- sum(rep(numberOfCases, fold) - foldLengths)

    correct <- c()
    correctLogistic <- c()
    probs <- list(length = numberOfCases)

    parents <- graph[[response]]
    variablesOfInterest <- c(response, parents)
    nlp <- sapply(parents, function(i) nlevels(data[,i]))
    nlr <- nlevels(data[,response])
    response_levels <- levels(data[, response])
    levels(data[, response]) <- paste("Level", levels(data[, response]),

                                      sep = ".")

    notresponse <- setdiff(seq_len(ncol(data)), response)

    # for priors??
    Nijk <- 1/prod(nlp)
    Nij <- sum(Nijk)

    numberIncorrect <- 0
    numberIncorrectLogistic <- 0
    numTests <- 0

    numberIncorrectByFold <- vector("numeric", fold)
    numberIncorrectLogisticByFold <- vector("numeric", fold)

    # loop over the N-folds of the cross validation
    for (foldNum in foldSeq){
      # set up the two dataset for this fold
      dataTest <- data[folds[[foldNum]], , drop = F]
      dataTraining <- data[-folds[[foldNum]], , drop = F]

      # number of data points in test data set
      numTestsFold <- nrow(dataTest)
      # ?
      numTests <- numTests + predictiveReplications * numTestsFold
      numberIncorrect <- 0
      numberIncorrectLogistic <- 0
      # cross tabulate the response against its parents
      # in the training data
      #
      # dataTrainingDF is a data frame
      # first columns are the parents

      createDataFrameTable <- function(data, response, conditioning){
        data <- table(data)
        data <- as.data.frame(data)

        data_molten <- melt(data, id.vars = names(data)[-length(names(data))])

        formula_response <- paste("`", names(data)[response], "`",
                                  collapse = "+", sep = "")
        formula_conditioning <- paste("`", names(data)[conditioning],"`",
                                      collapse = " + ", sep = "")

        formula_reshape <- eval(as.formula(
          paste(
            formula_conditioning,
            "~",
            formula_response
          )),

          data_molten
        )

        data_cast <- cast(
          data_molten,

          formula_reshape,
          sum,
          margins = "grand_col"
        )

        colnames(data_cast)[colnames(data_cast) == "(all)"] <- "total"

        data_cast
      }

      dataTrainingCast <- createDataFrameTable(dataTraining, response, parents)

      # an unfriendly version of transform.data.frame
      # some of this inspired by transform.data.frame
      mytransform <- function(prop_form, df){
        with(
        df,

        do.call(
          "data.frame",

          c(list(df), prop_form)
        )
      )
      }

      toparse <- paste("`Level.", response_levels, "` + Nijk", sep = "")
      posterior_formula <- as.list(parse(text = toparse))
      posterior_formula_names <- paste("Post.", response_levels, sep = "")
      names(posterior_formula) <- posterior_formula_names
      dataTrainingCast <- mytransform(posterior_formula, dataTrainingCast)

      ### Posterior totals
      toparse <- paste("`Post.", response_levels, "`", collapse = "+", sep = "")
      posterior_totals <- as.list(parse(text = toparse))
      posterior_totals_names <- "Post.total"
      names(posterior_totals) <- posterior_totals_names
      dataTrainingCast <- mytransform(posterior_totals, dataTrainingCast)

      ## Posterior probabilities
      toparse <- paste("`Post.", response_levels, "`/`Post.total`", sep = "")
      posterior_probabilities <- as.list(parse(text = toparse))
      posterior_probabilities_names <- paste("Prop.", response_levels, sep = "")
      names(posterior_probabilities) <- posterior_probabilities_names
      dataTrainingCast <- mytransform(posterior_probabilities, dataTrainingCast)

      dataTestCast <- createDataFrameTable(dataTest, response, parents)

      # iterate over rows of dataTestCast
      for (rowNum in seq_len(nrow(dataTestCast))){
        N <- dataTestCast[rowNum, "total"]

        # if this configuration was actually in the test data
        if (N > 0){
          probs <- dataTrainingCast[rowNum, posterior_probabilities_names]
          out <- tabulate(sample.int(nlr, N * predictiveReplications,
                                     replace = T, prob = probs), nbins = nlr)
          numberIncorrect <- numberIncorrect + sum(abs(out -

                            predictiveReplications *

                            dataTestCast[rowNum, paste("Level.",

                            response_levels, sep = "")]))/2

          if (verbose){
            cat("trainingcast:\n")
            print(dataTrainingCast[rowNum, ])
            cat("\n Probs:\n")
            print(probs)
            cat("\n out")
            print(out)
            cat("\n NumInc\n", numberIncorrect, "\n")
          }
        }
      }

      # it appears this handled an edge case where there is only one level.
      # is this still needed??
      #
      #responseNames <- names(data)[-response]
      #for (item in seq_along(names(dataTraining)[-response])){
      #  if (nlevels(dataTraining[,names(dataTraining)[-response][item]]) == 1){
      #    responseNames <- responseNames[-item]
      #  }
      #}

      makeQuotedFormulaString <- function(names){
        paste("`", names, "`", collapse = "+", sep = "")
      }

      form <- eval(
        paste(
          makeQuotedFormulaString(names(data)[response]),

          " ~ ",

          makeQuotedFormulaString(names(data)[-response])
        ),
        dataTraining
      )

      if (binomial){
        #lr <- glm(form, data = dataTraining, trace = F, family = "binomial")
        #numberIncorrectLogistic <-

        #  numberIncorrectLogistic +

        #  sum(rep(
        #    ifelse(
        #      predict(lr, dataTest[,-response, drop = F]) > 0,

        #      1,

        #      0
        #    ) != dataTest[, response],

        #    predictiveReplications
        #  ))
      }
      else {
        xyz <- as.data.frame(table(dataTraining))
        mult <- multinom(form, data = xyz, weights = Freq, trace = verbose)
        #browser()
        #numberIncorrectLogistic <- numberIncorrectLogistic +

        # sum(rep(predict(mult, dataTest[, - response]) !=

        # dataTest[, response], predictiveReplications))

        dataTestMult <- createDataFrameTable(dataTest, response, notresponse)

        for (rowNum in seq_len(nrow(dataTestMult))){
          N <- dataTestMult[rowNum, "total"]

          # if this configuration was actually in the test data
          if (N > 0){
            probs <- predict(mult,
                             dataTestMult[rowNum, names(data)[-response]],
                             type = "probs")
            # irritatingly predict.nnet doesn't
            # give a vector of probs for nlev(resp) = 2
            if (nlr == 2){
              probs <- c(1 - probs, probs)
            }
            #browser()
            out <- tabulate(sample.int(nlr, N * predictiveReplications,

                            replace = T, prob = probs), nbins = nlr)
            numberIncorrectLogistic <- numberIncorrectLogistic + sum(abs(out -

                                       predictiveReplications *

                                      dataTestMult[rowNum, paste("Level.",

                                      response_levels, sep = "")]))/2
          }
        }
      }

      if (verbose){
        cat(
          "BVS: After ", foldNum,

          " folds ", (numTests - numberIncorrect)/numTests,

          " were correctly predicted.\n", sep = "")
        cat(
          "Logisitic: After ", foldNum,

          " folds ", (numTests - numberIncorrectLogistic)/numTests,

          " were correctly predicted.\n", sep = "")
      }

      numberIncorrectLogisticByFold[foldNum] <- numberIncorrectLogistic
      numberIncorrectByFold[foldNum] <- numberIncorrect
    }
    if (verbose){
      cat("BVS:",
          sum(numberIncorrectByFold)/numTests,
          ". MultnomLR:",
          sum(numberIncorrectLogisticByFold)/numTests, "\n")
    }

    list(
      overall = c(bn = sum(numberIncorrectByFold),
                  lr = sum(numberIncorrectLogisticByFold)),
      numberIncorrectByFold = numberIncorrectByFold,
      numberIncorrectLogisticByFold = numberIncorrectLogisticByFold
    )
  } else {
    cat("nnet and reshape required")
  }
}
