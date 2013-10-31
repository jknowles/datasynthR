\name{rnegbinomcor}
\alias{rnegbinomcor}
\title{Generate a vector from a negative binomial distribution correlated with another distribution}
\usage{
  rnegbinomcor(x, rho, scalar = NULL)
}
\arguments{
  \item{x}{variable to draw from}

  \item{rho}{correlation coefficient between x and result
  of function}

  \item{scalar}{a scale factor for the negative binomial
  draws}
}
\value{
  a vector of the same length as x drawn from a negative
  binomial distribution correlated with x at the level of
  rho
}
\description{
  Allow user to draw from a negative binomial distribution
  correlated with a user specified vector
}
\details{
  Rough estimate
}
\author{
  Jared E. Knowles
}

