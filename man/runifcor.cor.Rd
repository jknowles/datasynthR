\name{runifcor.cor}
\alias{runifcor.cor}
\title{Generate a vector from a uniform distribution correlated with another distribution}
\usage{
  runifcor.cor(x, rho)
}
\arguments{
  \item{x}{variable to draw from}

  \item{rho}{correlation coefficient between x and result
  of function}
}
\value{
  a vector of the same length as x drawn from a normal
  distribution correlated with x at the level of rho
}
\description{
  Allow user to draw from a uniform distribution correlated
  with a user specified vector
}
\details{
  Rough estimate
}
\note{
  runif only works coming from uniform data
}
\author{
  Jared E. Knowles
}
\references{
  modified from Eric Neuwirth
  \url{http://r.789695.n4.nabble.com/Generating-uniformly-distributed-correlated-data-td3314905.html}
}

