\name{rbinomcor}
\alias{rbinomcor}
\title{Generate a vector from a binomial distribution correlated with another distribution}
\usage{
  rbinomcor(x, rho, scale = NULL)
}
\arguments{
  \item{x}{variable to draw from}

  \item{rho}{correlation coefficient between x and result
  of function}

  \item{scale}{a scale factor for the binomial draws}
}
\value{
  a vector of the same length as x drawn from a binomial
  distribution correlated with x at the level of rho
}
\description{
  Allow user to draw from a binomial distribution
  correlated with a user specified vector
}
\details{
  Rough estimate
}
\author{
  Jared E. Knowles
}

