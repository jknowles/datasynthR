##'
##'

# Build a correlated matrix

genNumeric <- function(n, k, names, rho, pattern){
  if(missing(pattern)){
  cov <- array(runif(n*k, -2, 2), dim=c(n, k))
  rnormcor <- function(x,rho) rnorm(1,rho*x,sqrt(1-rho^2))
  for(i in 2:k){
    cov[, i] <- sapply(cov[, i-1], rnormcor, rho=rho)
  }
  return(cov)
  } else if(!missing(pattern)){
    message("Code not written yet...")
  }
}

# Pass different distributional forms:

pattern <- list(dist = c("norm", "binom", "chisq", "pois", "unif", 
                         "weibull", "gamma"), 
                center = c(100, 0, 2, 4, 9, 2, 12), 
                spread = c())


# Generate glm


# evaluate name?
genBinomiaDV <- function(df, form, errors, intercept){
  # if form = NULL
  exp <- paste0(form$coefs, "*","gdf$", form$vars, collapse="+")
  form$exp <- parse(text=exp)
  z <- intercept + eval(form$exp) + rnorm(dim(gdf)[1])
  pr <- 1/(1+exp(-z))
  y = rbinom(dim(gdf)[1], 1, pr)
  return(y)
}


genFactor <- function(n, k, names, nlevel, rho, ...){
  cov <- genNumeric(n, k, names, rho)
  cov <- as.data.frame(cov)
  for(i in 1:ncol(cov)){
    cov[, i] <- cut(cov[, i], breaks=nlevel)
    cov[, i] <- factor(cov[, i], labels=sample(letters, nlevel))
  }
  return(as.data.frame(cov))
}
  
