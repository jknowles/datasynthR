################################################################################
# Test missingness
################################################################################

context("vector NA search")

test1 <- c(rnorm(100), rep(NA, 28))
test2 <- c(rnorm(1000), rep(NA, 7))
test3 <- c(rep(NA, 400), rnorm(20))

test_that("Correct randomness is detected", {
  expect_equivalent(vecNAsearch(test1), 28/128)
  expect_equivalent(vecNAsearch(test2), 7/1007)
  expect_equivalent(vecNAsearch(test3), 400/420)
})

test_that("Function is fast!", {
  expect_that(vecNAsearch(test1), takes_less_than(.002))
  expect_that(vecNAsearch(test2), takes_less_than(.002))
  expect_that(vecNAsearch(test3), takes_less_than(.002))
  
})



context("Vector MCAR")

P <- 0.1
N <- 1000
x <- rnorm(N)
y <- MCARx(x, P)
a <- length(y[is.na(y)])
t <- .05
e <- N * P
a2 <- abs(e - a) / N

test_that("Correct random chance of missingness is created", {
  expect_that(a2, is_less_than(t))
})

P <- 0.456
N <- 100000
x <- rnorm(N)
y <- MCARx(x, P)
a <- length(y[is.na(y)])
t <- .05
e <- N * P
a2 <- abs(e - a) / N

test_that("Function is fast and accurate", {
  expect_that(MCARx(x,P), takes_less_than(.02))
  expect_that(a2, is_less_than(t))
})


context("Matrix DF")


N <- 8000
K <- 10
P <- 0.1
RHO1 <- 0.45
covmat <- genNumeric(N, K, rho=RHO1)

covmat2 <- MCAR.df(covmat, P)

a <- apply(covmat2, 2, vecNAsearch) * N
t <- .05
e <- K * N * P
a2 <- abs(sum(a) - e) / N

test_that("Function is fast and accurate", {
  expect_that(MCAR.df(covmat,P), takes_less_than(.05))
  expect_that(a2, is_less_than(t))
})

N <- 98756
K <- 25
P <- 0.43278456
RHO2 <- -0.24
covmat <- genNumeric(N, K, rho=RHO2)

covmat2 <- MCAR.df(covmat, P)

a <- apply(covmat2, 2, vecNAsearch) * N
t <- .05
e <- K * N * P
a2 <- abs(sum(a) - e) / N

test_that("Function is fast and accurate", {
  expect_that(MCAR.df(covmat,P), takes_less_than(1))
  expect_that(a2, is_less_than(t))
})

context("dimNA function")


N <- 9874
K <- 25
P <- 0.43278456
RHO2 <- -0.24
covmat <- genNumeric(N, K, rho=RHO2)

covmat2 <- MCAR.df(covmat, P)

MR <- dimNA(covmat2)


#a <- apply(covmat2, 2, vecNAsearch) * N
t <- .005
e <- K * N * P
#a2 <- abs(sum(a) - e) / N

test_that("Function returns the correct dimensions answers", {
  expect_equivalent(MR$TotalCells, N*K)
  expect_that(abs(MR$TotalMissing - e) / (N*K), is_less_than(t))
  expect_that(abs(MR$TotalProportionMissing - P), is_less_than(t))
})

test_that("Function is fast!", {
  expect_that(dimNA(covmat2), takes_less_than(.1))
})

context("MAR data")


N <- 5000
K <- 25
P <- 0.43278456
RHO2 <- -0.24
covmat <- genNumeric(N, K, rho=RHO2)


covmat2 <- MAR(covmat, 1)

