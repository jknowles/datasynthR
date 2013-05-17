################################################################################
# Random variables with correlations
#
#################################################################################

rnormcor <- function(x,rho) rnorm(1, rho*x, sqrt(1-rho^2))

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


runifcor.cor <- function(x, rho){
  hw <- function(r){ 
    tmp <- (3-sqrt(1+8*abs(r)))/4 
    return(tmp * sign(r))
  } 
  y <- (x*sign(rho) + runif(1000,-hw(abs(rho)),hw(rho*sign(rho)))) %% 1 
  return(y)
}

# only works coming from uniform data



"unif", 
"weibull", "gamma"