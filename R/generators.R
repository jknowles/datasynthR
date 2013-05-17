##'
##'

# Build a correlated matrix

genNumeric <- function(n, k, names, rho, pattern){
  if(missing(pattern)){
  cov <- array(runif(n*k, -2, 2), dim=c(n, k))
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
                spread = c(3, 5, 2, 9, 4, 7, 3))


cov <- matrix(nrow = N, ncol = length(pattern$dist))

type <- match.arg(pattern$dist[1], c("norm", "binom", "chisq", "pois", "unif", 
                             "weibull", "gamma"))

cov[, 1] <- switch(type, 
                   norm = rnorm(cov[, 1], mean=pattern$center[1]),
                   binom = rbinom(cov[, 1], size=2*pattern$center[1], prob=0.5),
                   chisq = rchisq(cov[, 1], df=pattern$center[1]), 
                   pois = rpois(cov[, 1], lambda=pattern$center[1]), 
                   unif = runif(cov[, 1], min=pattern$center[1] - pattern$spread[1], 
                                max= pattern$center[1] + pattern$spread[1]),
                   weibull= rweibull(cov[, 1], scale=pattern$center[1], shape=pattern$spread[1]), 
                   gamma = rgamma(cov[, 1]), shape=pattern$center[1])
for(i in 2:ncol(cov)){
  type <- match.arg(pattern$dist[i], c("norm", "binom", "chisq", "pois", "unif", 
                                       "weibull", "gamma"))
  cov[, i]  <- switch(type, 
                      norm = rnorm(cov[, i], mean=pattern$center[i]),
                      binom = rbinom(cov[, i], size=2*pattern$center[i], prob=0.5),
                      chisq = rchisq(cov[, i], df=pattern$center[i]), 
                      pois = rpois(cov[, i], lambda=pattern$center[i]), 
                      unif = runif(cov[, i], min=pattern$center[i] - pattern$spread[i], 
                                   max= pattern$center[i] + pattern$spread[i]),
                      weibull= rweibull(cov[, i], scale=pattern$center[i], shape=pattern$spread[i]), 
                      gamma = rgamma(cov[, i], shape=pattern$spread[i]))
}


 
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
  
