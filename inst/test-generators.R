################################################################################
# Test generators
################################################################################

context("Simple Random Correlation Matrices")

# gdf <- genNumeric(8000, 50, rho=0.45)
# gdf <- as.data.frame(gdf)
# test <- genNumeric(1000, pattern=pattern)


N <- 8000
K <- 10
RHO1 <- 0.45
covmat <- genNumeric(N, K, rho=RHO1)

test_that("Dimensions are correct", {
  expect_equal(nrow(covmat), N)
  expect_equal(ncol(covmat), K)
  expect_is(covmat, "matrix")
  
})


covE <- rep(NA, K-1)

for(i in 2:ncol(covmat)){
  covE[i-1] <- cor(covmat)[i-1, i]
}

tol <- RHO1 * K * .005
err <- sum(abs(covE - RHO1)) / K

test_that("Covariances are correct", {
  expect_that(err, is_less_than(tol))
})






pattern <- list(dist = c("norm", "binom", "chisq", "pois", "unif", 
                         "weibull", "gamma"), 
                center = c(100, 0, 2, 4, 9, 2, 12), 
                spread = c(3, 5, 2, 9, 4, 7, 3), 
                rho = c(-.05, -.4, 0.3, 0.9, .03, -.6, -.2))


N <- 1000

cov <- matrix(nrow = 1000, ncol = length(pattern$dist))
