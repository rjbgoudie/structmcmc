\name{options.grid}
\alias{options.grid}
\title{input a list x.}
\usage{
  options.grid(x, maxIndegree, required = NULL)
}
\arguments{
  \item{x}{A list}

  \item{maxIndegree}{Maximum indegree}

  \item{required}{...}
}
\value{
  A list of options
}
\description{
  return a list that includes all options, including those
  of varying sizes
}
\details{
  eg x = list(c(1,2), c(2, 3)) out = list(list(c(1), c(2)),
  list(c(2), c(2)), list(c(1, 2), c(2)), list(c(1), c(3)),
  list(c(2), c(3)), list(c(1, 2), c(3)), list(c(1), c(2,
  3)), list(c(2), c(2, 3)), list(c(1, 2), c(2, 3)))
}

