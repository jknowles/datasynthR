\name{genNumeric}
\alias{genNumeric}
\title{Generate a dataframe of numeric variables with known distributions and relationships}
\usage{
  genNumeric(n, k, rho, seed, pattern)
}
\arguments{
  \item{n}{Number of rows of data to be generated}

  \item{k}{Number of columns to be generated}

  \item{rho}{Correlation coefficient between pairs of
  variables}

  \item{seed}{A vector of numerics length n to be used to
  generate correlations for other variabes from}

  \item{pattern}{List of attributes for columns of data in
  the data frame created}
}
\value{
  An R data frame of n rows and k columns with
  distributions specified in \code{\link{pattern}}. If
  \code{\link{pattern}} is not specified then variables are
  normally distributed with sequential bivariate
  correlations equal to rho.
}
\description{
  Quickly generate random numeric data with known
  properties
}
\details{
  \code{\link{pattern}} allows the user to specify a list
  with three elements: dist, rho, and name. Each element
  should be length k. Dist currently supports the options
  for norm (normal), binom (binomial), chisq (Chi-squared),
  pois (poisson), unif (uniform), weibull (Weibull), and
  gamma (gamma) distributions. Rho should be a numeric
  between -1 and 1 representing the correlation coefficient
  between that variable and the first column of the data
  frame. Names should be characters corresponding to the
  names of the columns in the resulting dataframe.
}
\note{
  For low n the value of rho will vary more greatly from
  the desired value.

  Currently, depending on the distribution the correlation
  is being built against, rho should be the correct sign,
  but it will not always result in the correct magnitude.
}
\examples{
dat1 <- genNumeric(1000, 3, rho=0.3)
cor(dat1[, 1], dat1[, 2])
cor(dat1[, 2], dat1[, 3])
# Specify a pattern
struc <- list(dist=c("norm", "pois", "unif"), rho=c(0.2, -.5, .5),
names=c("super", "cool", "data"))
dat2 <- genNumeric(1000, pattern=struc)
cor(dat2[, 1], dat2[, 2])
cor(dat2[, 1], dat2[, 3])
}
\author{
  Jared E. Knowles
}

