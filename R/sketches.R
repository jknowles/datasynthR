
#################
#
#
#
#

set.seed(214521)

rho <- 0.8

rnormcor <- function(x,rho) rnorm(1, rho*x, sqrt(1-rho^2))

x <- rnorm(10000) 

y <- sapply(x, rnormcor, rho=rho)

cor(x,y)

rnormcorV <- function(x, rho){
#   powerNeg <- function(x, p) (abs(x)^(1/p)) * sign(x)
#   dist <- rho - 0.5
#   steps <- qpois(dist, lambda=6)
  y <- sapply(x, rnormcor, rho=rho)
  return(y)
}

y2 <- rnormcorV(x, rho)

cor(x,y2)

y3 <- rnormcorV(x, 0.8)

cor(x, y3)

#http://www.r-bloggers.com/simulating-data-following-a-given-covariance-structure/

library(lattice) # for splom
library(car)     # for vif

# number of observations to simulate
nobs = 100

# Using a correlation matrix (let' assume that all variables
# have unit variance
M = matrix(c(1, 0.7, 0.7, 0.5,
             0.7, 1, 0.95, 0.3,
             0.7, 0.95, 1, 0.3,
             0.5, 0.3, 0.3, 1), nrow=4, ncol=4)

# Cholesky decomposition
L = chol(M)
nvars = dim(L)[1]

# R chol function produces an upper triangular version of L
# so we have to transpose it.
# Just to be sure we can have a look at t(L) and the
# product of the Cholesky decomposition by itself

t(L)

t(L) %*% L


# Random variables that follow an M correlation matrix
r = t(L) %*% matrix(rnorm(nvars*nobs), nrow=nvars, ncol=nobs)
r = t(r)

rdata = as.data.frame(r)
names(rdata) = c('resp', 'pred1', 'pred2', 'pred3')
