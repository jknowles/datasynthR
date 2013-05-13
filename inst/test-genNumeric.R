library(testthat)


gdf <- genNumeric(2000, 50, rho=0.05)

gdf <- as.data.frame(gdf)


terms <- list(coefs = rnorm(5), vars = sample(names(gdf), 5)) 

y <- genBinomiaDV(gdf, form=terms, intercept=3)