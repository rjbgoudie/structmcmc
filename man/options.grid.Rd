\name{options.grid}
\alias{options.grid}
\title{input a list x.}
\usage{options.grid(x)}

\description{
  input a list x.
}

\details{
  return a list that includes all options, including those
  of varying sizes

  eg x = list(c(1,2), c(2, 3)) out = list(list(c(1), c(2)),
  list(c(2), c(2)), list(c(1, 2), c(2)), list(c(1), c(3)),
  list(c(2), c(3)), list(c(1, 2), c(3)), list(c(1), c(2,
  3)), list(c(2), c(2, 3)), list(c(1, 2), c(2, 3)))
}
\value{A list of options}
\arguments{
  \item{x}{A list}
}
