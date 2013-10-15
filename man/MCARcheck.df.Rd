\name{MCARcheck.df}
\alias{MCARcheck.df}
\title{Test correlation of missingness for all pairwise combinations in a dataframe}
\usage{
  MCARcheck.df(df)
}
\arguments{
  \item{df}{a dataframe}
}
\value{
  a list with three matrices containing the pairwise Gamma
  statistic for all columns in the dataframe, the standard
  error of the Gamma statistic, and the significance test
  of the Gamma statistic
}
\description{
  Return the GK test for all pairs of columns in a
  dataframe on missing values
}
\details{
  This function forces each pair of columns to factors of
  missing and non-missing elements and then performs a GK
  test on those factors to determine if they are
  correlated, to what degree, and if this is statistically
  significant.
}
\note{
  Yadda yadda yadda
}
\author{
  Jared E. Knowles
}

