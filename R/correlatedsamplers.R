################################################################################
# Random variables with correlations
#
#################################################################################

rnormcor <- function(x,rho) rnorm(1,rho*x,sqrt(1-rho^2))

rchisqcor <- function(x, rho){
  sign(rho)*sign(x)*sum(rep(sapply(x, rnormcor, rho=rho), times=length(x))^2)
}

rpoiscor <- function(x, rho, scale){
  tmp <- sapply(x, rnormcor, rho=rho)
  y <- rpois(length(x), lambda=mean(tmp))
  return(y)
} 

a <- rnorm(1000)
tmp <- sapply(a, rnormcor, rho=0.7)
y <- sapply(tmp, rpois, lambda=3)

y <- sapply(tmp, rpois, lambda=(tmp+.01)^2)



b <- rpoiscor(a, rho=.1)

cor(a,b)


"unif", 
"weibull", "gamma"