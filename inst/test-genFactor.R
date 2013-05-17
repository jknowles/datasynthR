# test factors


context("Is a data.frame of factors created properly?")

N <- 1000
K <- 5
NLEVEL <- 5
RHO <- 0.80001
b <- genFactor(N, K, nlevel=NLEVEL, rho=RHO)

test_that("genFactor outputs data in expected size and shape", {
  expect_is(b, "data.frame")
  expect_equal(nrow(b), N)
  expect_equal(ncol(b), K)
  expect_equal(length(table(b[1])), NLEVEL)
})

z <- chisq.test(b[,1], b[,2], simulate.p.value=TRUE)
table(b[,1], b[,2])

context("Are the variables in the data.frame properly related?")

test_that("RHO is correct", {
  
  
  
})

context("Speed tests")

test_that("Function is not slow", {
  
})