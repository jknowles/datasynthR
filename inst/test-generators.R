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

test_that("Correlations are correct", {
  expect_that(err, is_less_than(tol))
})



context("User Controlled Random Correlation Matrices")


struc <- list(dist = c("norm", "binom", "chisq", "pois", "unif", 
                         "weibull", "gamma"), 
                rho = c(-.05, -.4, 0.3, 0.9, .03, -.6, -.2))


N <- 5000

covmat <- genNumeric(N, pattern=struc)

test_that("Dimensions are correct", {
  expect_equal(nrow(covmat), N)
  expect_equal(ncol(covmat), length(struc$dist))
  expect_is(covmat, "matrix")
  
})

err <- sum(abs(cor(covmat)[2:ncol(covmat),1] - struc$rho[2:length(struc$rho)])) / length(struc$rho)
tol <- sum(abs(struc$rho[2:length(struc$rho)]) *.05)

test_that("Correlations are correct", {
  expect_that(err, is_less_than(tol))
})

context("Speed")


N <- 98756
K <- 25
P <- 0.43278456
RHO2 <- -0.24
covmat <- genNumeric(N, K, rho=RHO1)

test_that("Function executes in reasonable length of time...", {
  expect_that(genNumeric(1000, 25, rho=0.3), takes_less_than(0.5))
  expect_that(genNumeric(500, 40, rho=0.3), takes_less_than(0.5))
  expect_that(genNumeric(87953, 8, rho=0.3), takes_less_than(0.5))
  expect_that(genNumeric(87953, 80, rho=0.3), takes_less_than(1))
})
