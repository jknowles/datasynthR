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



# evaluate name?
genBinomiaDV <- function(df, form, errors, intercept){
  # if form = NULL
  form$exp <- parse(text=paste0(form$coefs, "*","gdf$", form$vars))
  z <- intercept + eval(form$exp) + rnorm(dim(gdf)[1])
  pr <- 1/(1+exp(-z))
  y = rbinom(dim(gdf)[1], 1, pr)
  return(y)
}
