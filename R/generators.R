##' Generate a dataframe of numeric variables with known distributions and relationships
##'
##' Quickly generate random numeric data with known properties 
##' 
##' @param n Number of rows of data to be generated
##' @param k Number of columns to be generated
##' @param rho Correlation coefficient between pairs of variables
##' @param pattern List of attributes for columns of data in the data frame created
##' @return An R data frame of n rows and k columns with distributions specified in \code{\link{pattern}}.
##' If pattern is not specified then variables are normally distributed with sequential bivariate correlations
##' equal to rho. 
##' @details pattern allows the user to specify a list with three elements: dist, rho, and name. Each element should be
##' length k. Dist currently supports the options for norm (normal), binom (binomial), 
##' chisq (Chi-squared), pois (poisson), unif (uniform), weibull (Weibull), 
##' and gamma (gamma) distributions. Rho should be a numeric between -1 and 1 representing the 
##' correlation coefficient between that variable and the first column of the data frame. 
##' Names should be characters corresponding to the names of the columns in the resulting dataframe. 
##' @note For low n the value of rho will vary more greatly from the desired value.
##' @export
##' @author Jared E. Knowles
##' @note Currently, depending on the distribution the correlation is being built against, rho should be 
##' the correct sign, but it will not always result in the correct magnitude. 
##' @examples
##' dat1 <- genNumeric(1000, 3, rho=0.3)
##' cor(dat[, 1], dat[, 2])
##' cor(dat[, 2], dat[, 3])
##' # Specify a pattern
##' struc <- list(dist=c("norm", "pois", "unif"), rho=c(0.2, -.5, .5), 
##' names=c("super", "cool", "data"))
##' dat2 <- genNumeric(1000, pattern=struc)
##' cor(dat2[, 1], dat2[, 2])
##' cor(dat2[, 1], dat2[, 3])
genNumeric <- function(n, k, rho, pattern){
  if(missing(pattern)){
  cov <- array(runif(n*k, -2, 2), dim=c(n, k))
   for(i in 2:k){
    cov[, i] <- sapply(cov[, i-1], rnormcor, rho=rho)
  }
  return(cov)
  } else if(!missing(pattern)){
    cov <- matrix(nrow = n, ncol = length(pattern$dist))
    cov <- cbind(rnorm(nrow(cov)), cov)
    
    for(i in 2:ncol(cov)){
      type <- match.arg(pattern$dist[i-1], c("norm", "binom", "chisq", "pois", "unif", 
                                           "weibull", "gamma"))
      cov[, i]  <- switch(type, 
                          norm = rnormcorV(cov[, 1], rho=pattern$rho[i-1]),
                          binom = rbinomcor(cov[, 1], rho=pattern$rho[i-1]),
                          chisq = rchisqcor(cov[, 1], rho=pattern$rho[i-1]), 
                          pois = rpoiscor(cov[, 1], rho=pattern$rho[i-1]), 
                          unif = runifcor.cor(cov[, 1], rho=pattern$rho[i-1]),
                          weibull= rweibullcor(cov[, 1], rho=pattern$rho[i-1]), 
                          gamma = rgammacor(cov[, 1], rho=pattern$rho[i-1]))
    }
    if(!is.null(pattern$names)){
      cov <- as.data.frame(cov)
      names(cov) <- c("seed", pattern$names)
      return(cov)
    }
    return(as.data.frame(cov))
  }
}

##' Generate a binomial dependent variable
##'
# need some flexibility to norm the resulting coefficients to avoid all 1 and 0s
genBinomialDV <- function(df, form, errors, intercept){
  # if form = NULL
  exp <- paste0(form$coefs, "*","df$", form$vars, collapse="+")
  form$exp <- parse(text=exp)
  z <- intercept + eval(form$exp) + runif(dim(df)[1])
  pr <- 1/(1+exp(-z))
  y = rbinom(dim(df)[1], 1, pr)
  return(y)
}


##' Generate a dataframe of unordered factors with known associations
##' 
##' Quickly generate random dataframes of unordered factor variables with known associations
##' 
##' @param n Number of rows in the resulting dataframe
##' @param k Number of columns in the resulting dataframe
##' @param nlevel Number of levels in the factor variables created
##' @param rho Level of association among the variables created
##' 
##' @details Rho is used to generate associations between preceding variables in the dataframe. 
##' Element 1 and 2 are associated at the level of rho. Element 2 and 3 are also associated at the level of 
##' rho. All variables have the same number of levels -- nlevels -- and currently factor level names 
##' are randomly generated from \code{\link{letters}}.
##' @note For low n the value of rho will vary more greatly from the desired value.
##' @export
##' @author Jared E. Knowles
##' @examples
##' dat1 <- genFactor(1000, 12, nlevel=6, rho=0.4)
##' gammaGK(dat1[, 1], dat1[, 2])
##' gammaGK(dat1[, 2], dat1[, 3])
##' # Not close to Rho
##' gammaGK(dat1[, 1], dat1[, 3])
##' 
##' # low n deviates further from rho
##' dat2 <- genFactor(50, 10, nlevel=6, rho=0.2)
##' gammaGK(dat2[, 1], dat2[, 2])
genFactor <- function(n, k, nlevel, rho, ...){
  cov <- genNumeric(n, k, rho)
  cov <- as.data.frame(cov)
  for(i in 1:ncol(cov)){
    cov[, i] <- cut(cov[, i], breaks=nlevel)
    cov[, i] <- factor(cov[, i], labels=sample(letters, nlevel))
  }
  return(as.data.frame(cov))
}
  
