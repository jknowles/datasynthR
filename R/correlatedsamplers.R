################################################################################
# Random variables with correlations
#
#################################################################################

# consider adding shape modifiers 
# Noise function as well

rnormcor <- function(x,rho) rnorm(1, rho*x, sqrt(1-rho^2))

rnormcorV <- function(x, rho){
  tmp <- sapply(x, rnormcor, rho=rho)
  return(tmp)
}

rchisqcor <- function(x, rho){
  sign(rho)*sign(x)*sum(rep(sapply(x, rnormcor, rho=rho), times=length(x))^2)
}

rpoiscor <- function(x, rho){
  tmp <- sapply(x, rnormcor, rho=rho)
  b <- round(tmp,digits=0)
  b <- b + abs(min(b))
  b <- b^3
  return(b)
} 

# runif only works coming from uniform data
# modified from Eric Neuwirth here http://r.789695.n4.nabble.com/Generating-uniformly-distributed-correlated-data-td3314905.html

runifcor.cor <- function(x, rho){
  hw <- function(r){ 
    tmp <- (3-sqrt(1+8*abs(r)))/4 
    return(tmp * sign(r))
  } 
  y <- (x*sign(rho) + runif(1000,-hw(abs(rho)),hw(rho*sign(rho)))) %% 1 
  return(y)
}




rweibullcor <- function(x, rho) {
  require(MASS)
  y <- sapply(x, rnormcor, rho=rho)
  y2 <- pnorm(y)
  fit <- fitdistr(y2, densfun="weibull")
  y3 <- sapply(y2, qweibull, shape=fit$estimate[[1]], scale=fit$estimate[[2]])
}



rgammacor <- function(x, rho){
  require(MASS)
  y <- sapply(x, rnormcor, rho=rho)
  y2 <- pnorm(y)
  fit <- fitdistr(y2, densfun="gamma")
  y3 <- sapply(y2, qgamma, shape=fit$estimate[[1]], rate=fit$estimate[[2]])
}


rbinomcor <- function(x, rho){
  require(MASS)
  y <- sapply(x, rnormcor, rho=rho)
  y2 <- pnorm(y)
  #pr <- 1/(1+exp(-y2))
  y <- sapply(y2, qbinom, size=1, prob=abs(rho))
}



