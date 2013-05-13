##'
##'

# Build a correlated matrix

genNumeric <- function(n, k, names, rho, pattern){
  cov <- array(runif(n*k, -2, 2), dim=c(n, k))
  rnormcor <- function(x,rho) rnorm(1,rho*x,sqrt(1-rho^2))
  for(i in 2:k){
    cov[, i] <- sapply(cov[, i-1], rnormcor, rho=runif(1, min=-1, max=1))
  }
  return(cov)
}


# Generate glm

form <- list(coefs = rnorm(5), vars = sample(names(gdf), 5)) 

# evaluate name?
genBinomiaDV <- function(df, form, errors, intercept){
  z <- intercept + 
}

# 
# z <- 3 + 2 * cov[, 1] + 3 * cov[, 2] - 8 * cov[, 3] +
#   15 * cov[, 4] - 2 * cov[, 5]^2 + 12 * cov[, 27] + 
#   .01 * cov[, 12] + .5 * cov[, 32] + - .0321 * cov[, 27] + rnorm(1000)
# 
# pr = 1/(1+exp(-z))         
# y = rbinom(1000,1,pr)

