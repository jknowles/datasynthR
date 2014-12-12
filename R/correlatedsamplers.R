################################################################################
# Random variables with correlations
#
#################################################################################

##' Generate a single variable from a distribution correlated with another distribution
##'
##' Allow user to draw from a random normal distribution correlated with a user specified distribution
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @return a single numeric value drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rnormcor <- function(x, rho) rnorm(1, rho*x, sqrt(1-rho^2))

##' Generate a single variable from a distribution correlated with another distribution
##'
##' Allow user to draw from a random log-normal distribution correlated with a user specified distribution
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @return a single numeric value drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rlnormcor <- function(x, rho, ...) rlnorm(1, rho*x, sqrt(1-rho^2))

##' Generate a vector from a normal distribution correlated with another distribution
##'
##' Allow user to draw from a random normal distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate, biased by known amount for now
##' @author Jared E. Knowles
##' @export 
##' @examples
##' x <- rnorm(1000, 1, 1)
##' y <- rnormcorV(x, 0.2)
##' cor(x,y) # very close to 0.2
##' mean(y) # close to 0
##' sd(y)   # close to 1
rnormcorV <- function(x, rho, ...){
#   powerNeg <- function(x, p) (abs(x)^(1/p)) * sign(x)
#   dist <- rho - 0.5
#   steps <- qpois(dist, lambda=6)
  y <- sapply(x, rnormcor, rho=rho)
#   for(i in 1:steps){
#   y <- sapply(y, rnormcor, rho=powerNeg(rho, i))
#   }
#   
#   require(MASS)
#   y <- sapply(x, rnormcor, rho=rho)
#   y <- sapply(y, rnormcor, rho=sqrtNeg(abs(rho)))
#   y <- sapply(y, rnormcor, rho=sqrtNeg(sqrtNeg(abs(rho))))
#   y <- sapply(y, rnormcor, rho=sqrtNeg(sqrtNeg(sqrtNeg(abs(rho)))))
#   y2 <- pnorm(y)
#   fit <- fitdistr(y, densfun="normal")
#   y3 <- sapply(y2, qnorm, mean=fit$estimate[[1]], sd=fit$estimate[[2]])
#   # clean up extreme values
#   y3[!is.finite(y3) & sign(y3) == 1] <- max(y3[is.finite(y3)], na.rm=T)
#   y3[!is.finite(y3) & sign(y3) == -1] <- min(y3[is.finite(y3)], na.rm=T)
#   y3 <- scale(y3, center=TRUE, scale=TRUE)
  return(y)
}

##' Generate a vector from a chi-square distribution correlated with another distribution
##'
##' Allow user to draw from a random chi-square distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @param na.rm a logical indicating whether to fit the distribution excluding missing values 
##' or to fail on missing values
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
##' @examples
##' x <- rnorm(1000, 1, 1)
##' y <- rchisqcor(x, 0.2)
##' cor(x,y) # very close to 0.2
##' mean(y) 
##' sd(y)   
rchisqcor <- function(x, rho, na.rm = FALSE) {
  y <- sapply(x, rlnormcor, rho=rho)
  y2 <- plnorm(y)
  if(anyNA(y)){
    if(na.rm == TRUE){
      p <- vecNAsearch(y)
      warnmsg <- paste(p, "proportion of vector is missing and was excluded", sep = " ")
      warning(warnmsg)
      y <- na.omit(y)
      y <- sanitizeInf(y)
    } else if(na.rm == FALSE){
      stop("Missing values found, consider using na.rm=TRUE")
    }
  }
  fit <- fitdistr(y, densfun="chi-squared", list(df=2, ncp=1), 
                  upper=c(max(y) - 1, max(y) +1), lower=c(0.01, 1))
  y3 <- sapply(y2, qchisq, df=fit$estimate[[1]], ncp=fit$estimate[[2]])
  return(y3)
}


##' Generate a vector from a Poisson distribution correlated with another distribution
##'
##' Allow user to draw from a Poisson distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @param na.rm a logical indicating whether to fit the distribution excluding missing values 
##' or to fail on missing values
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
##' @examples
##' x <- rnorm(1000, 1, 1)
##' y <- rpoiscor(x, 0.2)
##' cor(x,y) # very close to 0.2
##' mean(y) 
##' sd(y)   
rpoiscor <- function(x, rho, na.rm = FALSE, ...){
  # consider allowing a meanshift to occur at user request
  y <- sapply(x, rlnormcor, rho=rho)
  y2 <- plnorm(y)
  y <- round(y, digits=0)
  if(anyNA(y) == TRUE){
    if(na.rm == TRUE){
      p <- vecNAsearch(y)
      warnmsg <- paste(p, "proportion of vector is missing and was excluded", sep = " ")
      warning(warnmsg)
      y <- na.omit(y)
      y <- sanitizeInf(y)
    } else if(na.rm == FALSE){
      stop("Missing values found, consider using na.rm=TRUE")
    }
  }
  fit <- fitdistr(y, densfun="Poisson")
  y3 <- sapply(y2, qpois, lambda=fit$estimate[[1]])
  return(y3)
} 


##' Generate a vector from a uniform distribution correlated with another distribution
##'
##' Allow user to draw from a uniform distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
##' @note runif only works coming from uniform data
##' @references modified from Eric Neuwirth \url{http://r.789695.n4.nabble.com/Generating-uniformly-distributed-correlated-data-td3314905.html}
##' @examples
##' x <- runif(1000)
##' y <- runifcor.cor(x, 0.2)
##' cor(x,y) # very close to 0.2
##' mean(y) 
##' sd(y)   
runifcor.cor <- function(x, rho){
  hw <- function(r){ 
    tmp <- (3-sqrt(1+8*abs(r)))/4 
    return(tmp * sign(r))
  } 
  y <- (x*sign(rho) + runif(1000,-hw(abs(rho)),hw(rho*sign(rho)))) %% 1 
  return(y)
}



##' Generate a vector from a Weibull distribution correlated with another distribution
##'
##' Allow user to draw from a Weibull distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @param na.rm a logical indicating whether to fit the distribution excluding missing values 
##' or to fail on missing values
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
##' @examples
##' x <- rnorm(1000)
##' y <- rweibullcor(x, 0.2)
##' cor(x,y) # very close to 0.2
##' mean(y) 
##' sd(y)  
rweibullcor <- function(x, rho, na.rm=FALSE) {
  y <- sapply(x, rnormcor, rho=rho)
  y2 <- pnorm(y)
  if(anyNA(y2)){
    if(na.rm == TRUE){
      p <- vecNAsearch(y2)
      warnmsg <- paste(p, "proportion of vector is missing and was excluded", sep = " ")
      warning(warnmsg)
      y2 <- na.omit(y2)
      y2 <- sanitizeInf(y2)
    } else if(na.rm == FALSE){
      stop("Missing values found, consider using na.rm=TRUE")
    }
  }
  fit <- try(fitdistr(y2, densfun="weibull"))
  y3 <- sapply(y2, qweibull, shape=fit$estimate[[1]], scale=fit$estimate[[2]])
}


##' Generate a vector from a Gamma distribution correlated with another distribution
##'
##' Allow user to draw from a Gamma distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @param na.rm a logical indicating whether to fit the distribution excluding missing values 
##' or to fail on missing values
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
##' @examples
##' x <- rnorm(1000)
##' y <- rgammacor(x, 0.2)
##' cor(x,y) # very close to 0.2
##' mean(y) 
##' sd(y)  
rgammacor <- function(x, rho, na.rm = FALSE){
  y <- sapply(x, rnormcor, rho=rho)
  y2 <- pnorm(y)
  if(anyNA(y2)){
    if(na.rm == TRUE){
      p <- vecNAsearch(y2)
      warnmsg <- paste(p, "proportion of vector is missing and was excluded", sep = " ")
      warning(warnmsg)
      y2 <- na.omit(y2)
      y2 <- sanitizeInf(y2)
    } else if(na.rm == FALSE){
      stop("Missing values found, consider using na.rm=TRUE")
    }
  }
  fit <- fitdistr(y2, densfun="gamma")
  y3 <- sapply(y2, qgamma, shape=fit$estimate[[1]], rate=fit$estimate[[2]])
}

##' Generate a vector from a negative binomial distribution correlated with another distribution
##'
##' Allow user to draw from a negative binomial distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @param scale a scale factor for the negative binomial draws
##' @param na.rm a logical indicating whether to fit the distribution excluding missing values 
##' or to fail on missing values
##' @return a vector of the same length as x drawn from a negative binomial distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rnegbinomcor <- function(x, rho, scale=NULL, na.rm = FALSE){
  if(missing(scale)){
    scale <- 10
  }
  y <- sapply(x, rnormcor, rho=rho)
  y2 <- pnorm(y)
  if(anyNA(y2)){
    if(na.rm == TRUE){
      p <- vecNAsearch(y2)
      warnmsg <- paste(p, "proportion of vector is missing and was excluded", sep = " ")
      warning(warnmsg)
      y2 <- na.omit(y2)
      y2 <- sanitizeInf(y2)
    } else if(na.rm == FALSE){
      stop("Missing values found, consider using na.rm=TRUE")
    }
  }
  fit <- fitdistr(ceiling(y2*scale), densfun="negative binomial")
  y3 <- sapply(y2, qnbinom, size=fit$estimate[[1]], mu=fit$estimate[[2]])
}

##' Generate a vector from a binomial distribution correlated with another distribution
##'
##' Allow user to draw from a binomial distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @param scale a scale factor for the binomial draws
##' @param na.rm a logical indicating whether to fit the distribution excluding missing values 
##' or to fail on missing values
##' @return a vector of the same length as x drawn from a binomial distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rbinomcor <- function(x, rho, scale=NULL, na.rm = FALSE){
  if(missing(scale)){
    scale <- 10
  }
  y <- rnegbinomcor(x, rho, scale, na.rm = na.rm)
  y2 <- y
  y2[y2 <= mean(y)] <- 0
  y2[y2 >= mean(y)] <- 1
  return(y2)
}


sanitizeInf <- function(x){
  if(anyInf(x) == TRUE){
    if(length(x) == length(x[!is.finite(x)])){
      stop("No non-finite values of x reported")
    }
    x[!is.finite(x)] <- sample(x[is.finite(x)], length(x[!is.finite(x)]), replace = TRUE)
    return(x)
  } else{
    return(x)
  }
}


## BETA, t, and geom, and regular binomial 