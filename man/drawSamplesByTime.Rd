\name{drawSamplesByTime}
\alias{drawSamplesByTime}
\title{Draw samples from a MCMC sampler, by time...}
\usage{drawSamplesByTime(sampler, time, burnin, burninleft, thin, verbose,
    samplesIncrement=10000)}
\description{Draw samples from a MCMC sampler, by time}
\details{Draws samples from an MCMC sampler for a specified number of seconds.}
\arguments{\item{sampler}{A sampler object}
\item{time}{The number of seconds to spend drawing samples.}
\item{burnin}{The number of samples to discard from the beginning of
the sample (includes previous runs)}
\item{burninleft}{The number of samples still to be discarded.}
\item{thin}{The frequency with which samples should be kept. eg for
\code{thin = 3}, every third sample will be kept.}
\item{verbose}{A logical. Should a progress bar be displayed?}
\item{samplesIncrement}{The size by which the samples should be
incremented.}}
