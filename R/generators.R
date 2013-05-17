##'
##'


genNumeric <- function(n, k, names, rho, pattern){
  if(missing(pattern)){
  cov <- array(runif(n*k, -2, 2), dim=c(n, k))
   for(i in 2:k){
    cov[, i] <- sapply(cov[, i-1], rnormcor, rho=rho)
  }
  return(cov)
  } else if(!missing(pattern)){
    cov <- matrix(nrow = n, ncol = length(pattern$dist))
    cov[, 1] <- rnorm(nrow(cov))
    
    for(i in 2:ncol(cov)){
      type <- match.arg(pattern$dist[i], c("norm", "binom", "chisq", "pois", "unif", 
                                           "weibull", "gamma"))
      cov[, i]  <- switch(type, 
                          norm = rnormcor(cov[, 1], rho=pattern$rho[i]),
                          binom = rbinomcor(cov[, 1], rho=pattern$rho[i]),
                          chisq = rchisqcor(cov[, 1], rho=pattern$rho[i]), 
                          pois = rpoiscor(cov[, 1], rho=pattern$rho[i]), 
                          unif = runifcor.cor(cov[, 1], rho=pattern$rho[i]),
                          weibull= rweibullcor(cov[, 1], rho=pattern$rho[i]), 
                          gamma = rgammacor(cov[, 1], rho=pattern$rho[i]))
    }
  }
  return(cov)
}

 
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
  
