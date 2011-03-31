\name{whichEdgesTogglable}
\alias{whichEdgesTogglable}
\title{Find togglable edges...}
\usage{whichEdgesTogglable(routes, constraintT)}
\description{Find togglable edges}
\details{Finds "edge-locations" in the graph that can be added or removed without
introducing a cycle into the graph}
\value{A logical matrix of the same dimension as the supplied matrices, 
with entries indicating whether the corresponding edge can be added or 
removed without introducing a cycle.}
\arguments{\item{routes}{The routes matrix of the network}
\item{constraintT}{The transpose of a constraint matrix}}
