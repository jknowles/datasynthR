\name{rpoiscor}
\alias{rpoiscor}
\title{Generate a vector from a Poisson distribution correlated with another distribution}
\usage{
  rpoiscor(x, rho)
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
  Allow user to draw from a Poisson distribution correlated
  with a user specified vector
}
\details{
  Rough estimate
}
\author{
  Jared E. Knowles
}

