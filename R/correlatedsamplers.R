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
rnormcor <- function(x,rho) rnorm(1, rho*x, sqrt(1-rho^2))

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
rlnormcor <- function(x,rho) rlnorm(1, rho*x, sqrt(1-rho^2))

##' Generate a vector from a normal distribution correlated with another distribution
##'
##' Allow user to draw from a random normal distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rnormcorV <- function(x, rho){
  tmp <- sapply(x, rnormcor, rho=rho)
  return(tmp)
}

##' Generate a vector from a chi-square distribution correlated with another distribution
##'
##' Allow user to draw from a random chi-square distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rchisqcor <- function(x, rho) {
  require(MASS)
  y <- sapply(x, rlnormcor, rho=rho)
  y2 <- plnorm(y)
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
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rpoiscor <- function(x, rho){
  # consider allowing a meanshift to occur at user request
  y <- sapply(x, rlnormcor, rho=rho)
  y2 <- plnorm(y)
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
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rweibullcor <- function(x, rho) {
  require(MASS)
  y <- sapply(x, rnormcor, rho=rho)
  y2 <- pnorm(y)
  fit <- fitdistr(y2, densfun="weibull")
  y3 <- sapply(y2, qweibull, shape=fit$estimate[[1]], scale=fit$estimate[[2]])
}


##' Generate a vector from a Gamma distribution correlated with another distribution
##'
##' Allow user to draw from a Gamma distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rgammacor <- function(x, rho){
  require(MASS)
  y <- sapply(x, rnormcor, rho=rho)
  y2 <- pnorm(y)
  fit <- fitdistr(y2, densfun="gamma")
  y3 <- sapply(y2, qgamma, shape=fit$estimate[[1]], rate=fit$estimate[[2]])
}

##' Generate a vector from a binomial distribution correlated with another distribution
##'
##' Allow user to draw from a binomial distribution correlated with a user specified vector
##' 
##' @param x variable to draw from
##' @param rho correlation coefficient between x and result of function
##' @return a vector of the same length as x drawn from a normal distribution correlated with x at the level of rho
##' @details Rough estimate
##' @author Jared E. Knowles
##' @export 
rbinomcor <- function(x, rho){
  require(MASS)
  y <- sapply(x, rnormcor, rho=rho)
  y2 <- pnorm(y)
  pr <- 1/(1+exp(-y2))
  y <- sapply(pr, qbinom, size=1, prob=abs(rho))
}