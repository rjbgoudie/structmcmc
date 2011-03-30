# Part of the "structural" package, http://github.com/rbtgde/structural
# 
# This software is distributed under the GPL-3 license.  It is free,
# open source, and has the attribution requirements (GPL Section 7) in
#   http://github.com/rbtgde/structural
# 
# Note that it is required that attributions are retained with each function.
#
# Copyright 2008 Robert J. B. Goudie, University of Warwick

context("gRain")

test_that("Basic test", {
  if (require(gRain)){
    d <- data.frame(
      a = factor(c(1, rep(3,2), rep(1, 7))),
      b = factor(c(2, rep(1, 4), rep(2, 5))),
      c = factor(c(2, rep(2, 3), rep(1, 6)))
    )

    net <- bn(integer(0), integer(0), c(1,2))
    out <- ml(net, d)

    g <- as.grain(out, net, d)

    # look at 'out' to see the probability tables of 'a' and 'b'

    exp <- list(a = array(data     = c(0.8, 0.2),
                          dim      = 2,
                          dimnames = list(a = c("1", "3"))))
    expect_that(querygrain(g, nodes = c("a"), type = "marginal"),
                equals(exp))

    exp <- list(b = array(data     = c(0.4, 0.6),
                          dim      = 2,
                          dimnames = list(b = c("1", "2"))))
    expect_that(querygrain(g, nodes = c("b"), type = "marginal"),
              equals(exp))

    # 0.62 = 0.8 * 0.4 * 0.5 +
    #        0.8 * 0.6 * 0.833333 +
    #        0.2 * 0.4 * 0 +
    #        0.2 * 0.6 * 0.5
    exp <- list(c = array(data     = c(0.62, 0.38),
                          dim      = 2,
                          dimnames = list(c = c("1", "2"))))
    expect_that(querygrain(g, nodes = c("c"), type = "marginal"),
              equals(exp))

    # 2x2 table
    exp <- array(data     = c(0.32, 0.08, 0.48, 0.12),
                 dim      = c(2, 2),
                 dimnames = list(a = c("1", "3"), b = c("1", "2")))
    expect_that(querygrain(g, nodes = c("a", "b"), type = "joint"),
                equals(exp))

    # full joint table
    exp <- array(data     = c(0.16, 0, 0.4, 0.06, 0.16, 0.08, 0.08, 0.06),
                 dim      = c(2, 2, 2),
                 dimnames = list(a = c("1", "3"),
                                 b = c("1", "2"),
                                 c = c("1", "2")))
    expect_that(querygrain(g, nodes = c("a", "b", "c"), type = "joint"),
                equals(exp))

    # conditioning
    g2 <- setFinding(g, nodes = c("a", "b"), states = c("1", "1"))

    exp <- list(c = array(data     = c(0.5, 0.5),
                          dim      = 2,
                          dimnames = list(c = c("1", "2"))))
    expect_that(querygrain(g2, nodes = c("c"), type = "marginal"),
                equals(exp))

    g2 <- setFinding(g, nodes = c("a", "b"), states = c("1", "2"))   
    exp <- list(c = array(data     = c(5/6, 1/6),
                          dim      = 2,
                          dimnames = list(c = c("1", "2"))))
    expect_that(querygrain(g2, nodes = c("c"), type = "marginal"),
                equals(exp))

    # reverse conditioning
    # from bayes rule
    #
    g2 <- setFinding(g, nodes = c("c"), states = c("1"))
    exp <- list(a = array(data     = c(0.9032258, 0.0967742),
                          dim      = 2,
                          dimnames = list(a = c("1", "3"))))
    expect_that(querygrain(g2, nodes = c("a"), type = "marginal"),
                equals(exp))

  } else {
    cat("gRain not installed")
  }
})

test_that("More advanced", {

  net <- structure(list(integer(0), integer(0), integer(0), integer(0),
      integer(0), integer(0), integer(0), c(3L, 9L), c(23L, 33L
      ), c(4L, 27L), c(33L, 35L), c(9L, 33L), c(1L, 33L), c(7L,
      32L), c(2L, 4L, 33L), c(22L, 37L), integer(0), c(3L, 33L,
      35L), c(1L, 3L, 18L), c(12L, 28L), c(4L, 6L, 33L), c(3L,
      18L), integer(0), c(16L, 25L), c(4L, 30L, 36L), c(8L, 27L
      ), 8L, c(1L, 21L), c(1L, 12L, 33L), c(3L, 35L), c(3L, 30L
      ), 30:31, c(1L, 3L, 6L), 24:25, c(4L, 17L, 33L), 4L, c(8L,
      26L), c(8L, 52L), c(9L, 37L), c(3L, 10L), c(18L, 33L, 52L
      ), c(12L, 49L), c(13L, 42L), 39L, c(15L, 41L), 16L, c(18L,
      33L, 37L), c(19L, 47L), c(20L, 52L), c(21L, 28L), c(1L, 29L,
      42L), c(26L, 37L)), class = c("bn", "parental"))

  dat <- structure(list(w1.BIO_SEX = structure(c(1L, 1L, 1L, 1L, 1L, 1L,
    1L, 2L, 2L, 2L), .Label = c("Male", "Female"), class = "factor"),
    w1.H1GI4 = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L), .Label = c("Not Hispanic/Spanish", "Hispanic/Spanish"
    ), class = "factor"), w1.H1GI6A = structure(c(1L, 1L, 2L,
    2L, 2L, 1L, 1L, 1L, 2L, 2L), .Label = c("Not white", "White"
    ), class = "factor"), w1.H1GI6B = structure(c(2L, 2L, 1L,
    1L, 1L, 2L, 2L, 2L, 1L, 1L),
    .Label = c("Not Black or African American",
    "Black or African American"), class = "factor"),
    w1.H1GI6C = structure(c(1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    .Label = c("Not American Indian or Native American",
    "American Indian or Native American"), class = "factor"),
    w1.H1GI6D = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L),
    .Label = c("Not Asian or Pacific Islander", "Asian or Pacific Islander"
    ), class = "factor"), w1.H1GI6E = structure(c(1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Not Other race",
    "Other race"), class = "factor"), w1.H1TO18 = structure(c(1L,
    1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Never",
    "Infrequently", "Often", "Almost daily"), class = "factor"),
    w1.H1ED2 = structure(c(1L, 1L, 1L, 3L, 1L, 1L, 1L, 1L, 2L,
    1L), .Label = c("Never", "A couple of times", "3-10 times",
    "11 or more"), class = "factor"), w1.H1ED21 = structure(c(2L,
    3L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L), .Label = c("Strongly disagree",
    "Neither strongly agree nor disagree", "Strongly agree"),
    class = "factor"),
    w1.H1DS5 = structure(c(3L, 2L, 3L, 3L, 2L, 2L, 1L, 1L, 1L,
    1L), .Label = c("Never", "1 or 2 times", "3 or 4 times",
    "5 or more times"), class = "factor"), w1.H1GH26 = structure(c(1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("No", "Yes"
    ), class = "factor"), w1.H1GH54 = structure(c(1L, 1L, 1L,
    1L, 2L, 2L, 1L, 1L, 1L, 1L), .Label = c("Minor", "Serious",
    "Extremely serious"), class = "factor"), w1.H1CO16D = structure(c(1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("No", "Yes"
    ), class = "factor"), w1.H1FV1 = structure(c(2L, 2L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Never", "Once",
    "More than once"), class = "factor"), w1.H1PF1 = structure(c(4L,
    4L, 4L, 3L, 4L, 3L, 4L, 4L, 4L, 4L),
    .Label = c("No resident mother figure",
    "Disagree", "Neither strongly agree or disagree", "Strongly agree"
    ), class = "factor"), w1.PC38 = structure(c(1L, 2L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L), .Label = c("No", "Yes"), class = "factor"),
    w1.H1ED7 = structure(c(2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 1L,
    1L), .Label = c("No", "Yes"), class = "factor"), w1.H1ED9 =
    structure(c(1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("No", "Yes"
    ), class = "factor"), w1.H1GH1 = structure(c(2L, 2L, 3L,
    3L, 3L, 3L, 2L, 3L, 2L, 3L), .Label = c("Poor", "Good/Fair",
    "Excellent/Very good"), class = "factor"), w1.H1NB2 = structure(c(2L,
    2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), .Label = c("False",
    "True"), class = "factor"), w1.PC34B = structure(c(4L, 3L,
    4L, 2L, 3L, 4L, 3L, 5L, 4L, 5L), .Label = c("Never", "Seldom",
    "Sometimes", "Often", "Always"), class = "factor"), w1.AGE =
    structure(c(1L,
    4L, 2L, 2L, 2L, 3L, 3L, 1L, 1L, 3L), .Label = c("11-12",
    "13-14", "15-16", "17-18", "19-21"), class = "factor"),
    w1.LIVEMOTHER = structure(c(2L,
    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L),
    .Label = c("Don't live with mother",
    "Live with mother"), class = "factor"),
    w1.LIVEFATHER = structure(c(2L,
    2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 2L),
    .Label = c("Don't live with father",
    "Live with father"), class = "factor"),
    w1.SMOKING = structure(c(1L,
    1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 1L), .Label = c("Non-smoker",
    "Former smoker", "Seven or fewer", "More than seven"),
    class = "factor"),
    w1.ALCOHOL = structure(c(1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L,
    1L, 1L), .Label = c("Never", "Once ever to once a month",
    "Twice a month to twice a week", "3 days a week to every day"
    ), class = "factor"), w1.EXERCISE = structure(c(3L, 2L, 2L,
    2L, 3L, 1L, 2L, 3L, 2L, 2L), .Label = c("No exercise", "Some exercise",
    "A lot of exercise"), class = "factor"), w1.CESD = structure(c(1L,
    2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Not depressive",
    "Depressive symptoms"), class = "factor"), w1.POVERTY = structure(c(2L,
    5L, 4L, 5L, 3L, 4L, 5L, 4L, 5L, 2L), .Label = c("5", "4",
    "3", "2", "1"), class = "factor"), w1.PARENT_ALC = structure(c(2L,
    2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L), .Label = c("Never",
    "Infrequent", "Regular", "More than 5 drinks"), class = "factor"),
    w1.HH_SMOKING = structure(c(2L, 3L, 3L, 2L, 1L, 1L, 2L, 1L,
    1L, 1L), .Label = c("No smokers", "Smoker at home", "Parent smokes"
    ), class = "factor"), w1.VIOL_VICTI = structure(c(2L, 2L,
    2L, 2L, 1L, 1L, 2L, 1L, 1L, 1L), .Label = c("Not victim of violence",
    "Victim of violence"), class = "factor"), w1.FAM_DEATH = structure(c(1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("No bereavement",
    "1 bereavement", "2+ bereavements"), class = "factor"),
    w1.ACAD_QUART = structure(c(2L,
    1L, 1L, 1L, 2L, 3L, 1L, 3L, 4L, 1L), .Label = c("Bottom",
    "Lower middle", "Upper middle", "Top"), class = "factor"),
    w1.PARNTS_REL = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 4L, 1L,
    3L, 1L), .Label = c("Parent relationship stable",
    "Argue but not talked of separating",
    "Talked of separating", "Talked of separating and argue a lot"
    ), class = "factor"), w1.DRUG_USER = structure(c(1L, 1L,
    1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Never",
    "Have used an illegal drug"
    ), class = "factor"), w2.H2TO22 = structure(c(1L, 1L, 1L,
    2L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Never", "Infrequently",
    "Often", "Almost daily"), class = "factor"),
    w2.H2ED2 = structure(c(1L,
    1L, 1L, 3L, 1L, 2L, 1L, 1L, 1L, 1L), .Label = c("Never",
    "A couple of times", "3-10 times", "11 or more"), class = "factor"),
    w2.H2ED17 = structure(c(1L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 2L,
    2L), .Label = c("Strongly disagree",
    "Neither strongly agree nor disagree",
    "Strongly agree"), class = "factor"), w2.H2FV16 = structure(c(3L,
    2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Never",
    "1 or 2 times", "3 or 4 times", "5 or more times"), class = "factor"),
    w2.H2GH28 = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L,
    1L), .Label = c("No", "Yes"), class = "factor"), w2.H2GH47 =
    structure(c(1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 1L), .Label = c("Minor",
    "Serious", "Extremely serious"), class = "factor"),
    w2.H2CO19D = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    .Label = c("No", "Yes"), class = "factor"),
    w2.H2FV1 = structure(c(3L, 2L, 1L,
    2L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Never", "Once",
    "More than once"), class = "factor"), w2.H2PF1 = structure(c(3L,
    3L, 2L, 3L, 3L, 4L, 3L, 4L, 4L, 4L),
    .Label = c("No resident mother figure",
    "Disagree", "Neither strongly agree or disagree", "Strongly agree"
    ), class = "factor"), w2.H2ED3 = structure(c(1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("No", "Yes"), class = "factor"),
    w2.H2ED5 = structure(c(2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L), .Label = c("No", "Yes"), class = "factor"),
    w2.H2GH1 = structure(c(2L,
    3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 2L), .Label = c("Poor", "Good/Fair",
    "Excellent/Very good"), class = "factor"), w2.H2NB2 =
    structure(c(2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L),
    .Label = c("False", "True"), class = "factor"),
    w2.CESD = structure(c(1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    .Label = c("Not depressive",
    "Depressive symptoms"), class = "factor"), w2.DRUG_USER =
    structure(c(1L,
    1L, 1L, 2L, 1L, 1L, 1L, 1L, 2L, 1L), .Label = c("Never",
    "Have used an illegal drug"), class = "factor")), .Names = c("w1.BIO_SEX",
    "w1.H1GI4", "w1.H1GI6A", "w1.H1GI6B", "w1.H1GI6C", "w1.H1GI6D",
    "w1.H1GI6E", "w1.H1TO18", "w1.H1ED2", "w1.H1ED21", "w1.H1DS5",
    "w1.H1GH26", "w1.H1GH54", "w1.H1CO16D", "w1.H1FV1", "w1.H1PF1",
    "w1.PC38", "w1.H1ED7", "w1.H1ED9", "w1.H1GH1", "w1.H1NB2", "w1.PC34B",
    "w1.AGE", "w1.LIVEMOTHER", "w1.LIVEFATHER", "w1.SMOKING", "w1.ALCOHOL",
    "w1.EXERCISE", "w1.CESD", "w1.POVERTY", "w1.PARENT_ALC", "w1.HH_SMOKING",
    "w1.VIOL_VICTI", "w1.FAM_DEATH", "w1.ACAD_QUART", "w1.PARNTS_REL",
    "w1.DRUG_USER", "w2.H2TO22", "w2.H2ED2", "w2.H2ED17", "w2.H2FV16",
    "w2.H2GH28", "w2.H2GH47", "w2.H2CO19D", "w2.H2FV1", "w2.H2PF1",
    "w2.H2ED3", "w2.H2ED5", "w2.H2GH1", "w2.H2NB2", "w2.CESD", "w2.DRUG_USER"
    ), row.names = c(NA, 10L), class = "data.frame")

  out <- bayes(net, dat, prior = "qi")

  g <- as.grain(out, net, dat)

  querygrain(g, nodes = c("w2.CESD"), type = "joint")

  g2 <- setFinding(g, nodes = c("w1.BIO_SEX"), states = c("Male"))
  querygrain(g2, nodes = c("w2.CESD"), type = "joint")

  g2 <- setFinding(g, nodes = c("w1.BIO_SEX"), states = c("Female"))
  querygrain(g2, nodes = c("w2.CESD"), type = "joint")

  z1 <- marginalGivenOthers(g, "w2.CESD", dat, FUN = function(x) x[2])

  z2 <- do.call("make.groups", z1)

  z2[, "exactlyWhich"] <- strtrim(rownames(z2), 20)
  dotplot(reorder(exactlyWhich, data) ~ data,
          data = z2[order(z2$data, decreasing = T)[1:20],])

  dotplot(reorder(exactlyWhich, data) ~ data,
          data = z2,
          scales = list(y = list(cex = 0.3)))

  z2[, "exactlyWhich"] <- ""
  dotplot(data ~ exactlyWhich | which,
          data = z2,
          scales = list(x = list(rot = 90)),
          horizontal = F)

  dotplot(data ~ which,
          data = z2,
          scales = list(x = list(rot = 90)),
          horizontal = F)

  z3 <- unlist(lapply(z1, function(node){
    diff(range(unlist(node)))
  }))
  z3df <- as.data.frame(z3)
  z3df[, "name"] <- rownames(z3df)

  dotplot(z3 ~ name,
          data = z3df,
          scales = list(x = list(rot = 90)))

})