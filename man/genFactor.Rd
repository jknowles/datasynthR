\name{genFactor}
\alias{genFactor}
\title{Generate a dataframe of unordered factors with known associations}
\usage{
  genFactor(n, k, nlevel, rho, ...)
}
\arguments{
  \item{n}{Number of rows in the resulting dataframe}

  \item{k}{Number of columns in the resulting dataframe}

  \item{nlevel}{Number of levels in the factor variables
  created}

  \item{rho}{Level of association among the variables
  created}
}
\description{
  Quickly generate random dataframes of unordered factor
  variables with known associations
}
\details{
  Rho is used to generate associations between preceding
  variables in the dataframe. Element 1 and 2 are
  associated at the level of rho. Element 2 and 3 are also
  associated at the level of rho. All variables have the
  same number of levels -- nlevels -- and currently factor
  level names are randomly generated from
  \code{\link{letters}}.
}
\note{
  For low n the value of rho will vary more greatly from
  the desired value.
}
\examples{
dat1 <- genFactor(1000, 12, nlevel=6, rho=0.4)
gammaGK(dat1[, 1], dat1[, 2])
gammaGK(dat1[, 2], dat1[, 3])
# Not close to Rho
gammaGK(dat1[, 1], dat1[, 3])

# low n deviates further from rho
dat2 <- genFactor(50, 10, nlevel=6, rho=0.2)
gammaGK(dat2[, 1], dat2[, 2])
}
\author{
  Jared E. Knowles
}

