##' Generate a dataframe of numeric variables with known distributions and relationships
##'
##' Quickly generate random numeric data with known properties 
##' 
##' @param n Number of rows of data to be generated
##' @param k Number of columns to be generated
##' @param rho Correlation coefficient between pairs of variables
##' @param pattern List of attributes for columns of data in the data frame created
##' @param seed A vector of numerics length n to be used to generate correlations for other variabes from
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
genNumeric <- function(n, k, rho, seed, pattern){
  if(missing(seed)){
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
  } else if(!missing(seed)){
    if(missing(pattern)){
      if(length(seed) != n) stop("Constant seed produces nonsense correlations.")
      cov <- array(runif(n*k, -2, 2), dim=c(n, k))
      cov <- cbind(seed, cov)
      for(i in 1:k+1){
        cov[, i] <- sapply(cov[, 1], rnormcor, rho=rho)
      }
      return(cov)
    } else if(!missing(pattern)){
      cov <- matrix(nrow = n, ncol = length(pattern$dist))
      #cov <- cbind(rnorm(nrow(cov)), cov)
      
      for(i in 1:ncol(cov)){
        type <- match.arg(pattern$dist[i], c("norm", "binom", "chisq", "pois", "unif", 
                                               "weibull", "gamma"))
        cov[, i]  <- switch(type, 
                            norm = rnormcorV(pattern$seed[,i], rho=pattern$rho[i]),
                            binom = rbinomcor(pattern$seed[,i], rho=pattern$rho[i]),
                            chisq = rchisqcor(pattern$seed[,i], rho=pattern$rho[i]), 
                            pois = rpoiscor(pattern$seed[,i], rho=pattern$rho[i]), 
                            unif = runifcor.cor(pattern$seed[,i], rho=pattern$rho[i]),
                            weibull= rweibullcor(pattern$seed[,i], rho=pattern$rho[i]), 
                            gamma = rgammacor(pattern$seed[,i], rho=pattern$rho[i]))
      }
      if(!is.null(pattern$names)){
        cov <- as.data.frame(cov)
        names(cov) <- c(pattern$names)
        return(cov)
      }
      return(as.data.frame(cov))
    }
    
  }
}

##' Generate a binomial dependent variable
##'
##' Allow the user to specify a formula for generating a binomial dependent variable
##' 
##' @param df dataframe to generate the dependent variable from
##' @param form A named list containing the variables in df to be used to calculate the binomial probability, and 
##' coefficients of those variables
##' @param errors The way the error structure of the DGP formula should be established
##' @param intercept An adjustment to the base probability
##' @return A binomial vector by a formula generated out of the elements of form
##' @details Coefficients needs to be long enough to incorporate the factor levels
##' @note Yadda yadda yadda
##' @export
##' @author Jared E. Knowles
##' @note Currently it can be easy for the user to build a formulae that results in all 0 or 1 results. 
##' Use intercept to adjust accordingly. Additionally, coefficient scaes don't make sense at the moment. 
##' Still need to add the ability to have confounders in place.
genBinomialDV <- function(df, form, errors, intercept, type = c("binary", "response")){
  if (missing(type)){
    type <- "binary"
  } else {
    type <- match.arg(type)
  }
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  exp <- paste0(form$coefs, "*","mod$", genFormula(df, form$vars)[-1], collapse=" + ")
  form$exp <- parse(text=exp)
  mm <- eval(parse(text=paste0("terms(~", paste0(form$vars, collapse=' + '),")")))
  mod <- as.data.frame(model.matrix(mm, df))
  mod$z <- intercept + eval(form$exp) + runif(dim(mod)[1])
  mod$pr <- 1/(1+exp(-mod$z))
  mod$pr <- range01(scale(mod$pr, center=TRUE))
  mod$y <- rbinom(dim(mod)[1], 1, mod$pr)
  if(type=="binary"){
    return(mod$y)
  } else if(type=="response"){
    return(mod$pr)
  }
}

##' Generate model matrix terms 
##'
##' Allow the user to determine the number of coefficients needed 
##' 
##' @param df dataframe to generate the dependent variable from
##' @param vars A vector of variable names in df to generate the formula terms from
##' @return A character vector representing the terms in a formula resulting from using the 
##' terms in vars
##' @details Help figure out
##' @note Yadda yadda yadda
##' @export
##' @author Jared E. Knowles
##' @note Convenience  function to help user understand what factor terms need to be expanded 
genFormula <- function(df, vars){
  z <- eval(parse(text=paste0("terms(~", paste0(vars, collapse=' + '),")")))
  test <- model.matrix(z, df)
  return(dimnames(test)[[2]])
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
genFactor <- function(n, k, nlevel, rho, seed, keepSeed, ...){
  if(missing(keepSeed)){
    keepSeed <- TRUE
  } else {keepSeed <- keepSeed}
  if(nlevel >= length(letters)){
    zz <- expand.grid(letters, LETTERS)
    smp <- unique(paste0(zz[,1], zz[,2])) #
  } else if(nlevel < length(letters)) {
    smp <- letters
  }
  if(missing(seed)){
  cov <- genNumeric(n, k, rho, ...)
  cov <- as.data.frame(cov)
  for(i in 1:ncol(cov)){
    cov[, i] <- cut(cov[, i], breaks=nlevel)
    draw <- length(unique(cov[, i]))
    cov[, i] <- factor(cov[, i], labels=sample(draw, nlevel))
  }
  if(keepSeed == TRUE){
    return(as.data.frame(cov))
  } else if(keepSeed == FALSE){
    cov <- cov[, 2:ncol(cov)]
    return(as.data.frame(cov))
  }
  }
  if(!missing(seed)){
    if(class(seed)=="numeric"){
    cov <- genNumeric(n, k, rho, seed)
    cov <- as.data.frame(cov)
    for(i in 1:ncol(cov)){
      cov[, i] <- cut(cov[, i], breaks=nlevel)
      draw <- length(unique(cov[, i]))
      cov[, i] <- factor(cov[, i], labels=sample(smp, draw))
    }
    } else if(class(seed)!="numeric"){
      seed.tmp <- as.numeric(as.factor(seed))
      cov <- genNumeric(n, k, rho, seed.tmp)
      cov <- as.data.frame(cov)
      cov[, 1] <- seed
      for(i in 2:ncol(cov)){
        cov[, i] <- cut(cov[, i], breaks=nlevel)
        draw <- length(unique(cov[, i]))
        cov[, i] <- factor(cov[, i], labels=sample(smp, draw))
      }
    }
    if(keepSeed == TRUE){
      return(as.data.frame(cov))
    } else if(keepSeed == FALSE){
      cov <- cov[, 2:ncol(cov)]
      return(as.data.frame(cov))
    }
  }
}
  # check errors for very large factor levels
