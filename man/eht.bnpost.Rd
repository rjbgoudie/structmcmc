\name{eht.bnpost}
\alias{eht.bnpost}
\title{Expected hitting time.}
\usage{\method{eht}{bnpost} (x, head=5, tp, ...)}
\description{Expected hitting time.}
\details{Computes the expected hitting times to the posterior modal graph from
the top 'head' other graphs.

See Norris (1998) Markov Chains. Cambridge. Theorem 1.3.5, p17.}
\value{A vector of the expected hitting times.}
\arguments{\item{x}{An object of class "bnpost".}
\item{head}{The expected hitting time from the top 'head' ranked graphs,
with respect to the posterior will be returned. All are computed,
so this setting has no effect on the run-time.}
\item{tp}{Optionally provide the transition probability matrix.}
\item{...}{Further arguments (unused)}}
