library(testthat)


gdf <- genNumeric(8000, 50, rho=0.45)

gdf <- as.data.frame(gdf)



# 
# 
# 
# terms <- list(coefs = runif(5, min=-3, max=3), vars = sample(names(gdf), 5)) 
# 
# y <- genBinomiaDV(gdf, form=terms, intercept=0)
# 
# 
# gdf <- cbind(y, gdf)
# 
# fullmod <- glm(y~., data=gdf, family="binomial")
# 
# form <- paste(terms$vars, collapse="+")
# form <- paste("y~", form)
# form <- as.formula(form)
# 
# 
# testmod <- glm(form, data=gdf, family="binomial")
# 
