################################################################################
# Test data generating process

set.seed(3240)

context("Binomial GLM")


N <- 8000
K <- 25
RHO1 <- 0.65
gdf <- genNumeric(N, K, rho=RHO1)
gdf <- as.data.frame(gdf)

nvars <- 3
modterms <- list(coefs = runif(nvars, min=-1, max=1), vars = sample(names(gdf), nvars)) 
 
y <- genBinomialDV(gdf, form=modterms, intercept=0)
gdf <- cbind(y, gdf)


test_that("Result are right class and dimensions", {
  expect_is(gdf, "data.frame")
  expect_equal(nrow(gdf), N)
  expect_equal(ncol(gdf), K+1)
})



fullmod <- glm(y~., data=gdf, family="binomial")
 
form <- paste(modterms$vars, collapse="+")
form <- paste("y~", form)
form <- as.formula(form)
testmod <- glm(form, data=gdf, family="binomial")


summary(testmod)
summary(fullmod)
anova(fullmod, testmod, test="LRT")
