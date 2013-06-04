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


test_that("Function produces data MCAR", {
  
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

set.seed(382)
seeds <- genNumeric(10000, 6, rho=0.1)
struc <- list(dist=c("norm", "norm", "unif", "pois", "pois", "gamma", 
                     "weibull"), 
              rho=c(0.7, 0.3, -0.5, 0.3, -0.8, 0.05, 0.7), 
              names=c("test1", "test2", "noise", "daysattended", 
                      "daysOUT", "bad", "bad2"), 
              seed = cbind(seeds[,1], seeds[,2], seeds[,3], seeds[, 4], seeds[, 5], 
                           seeds[, 6], seeds[,1]))
dat <- genNumeric(10000, pattern=struc)
dat1 <- genFactor(10000, 3, nlevel=3, rho=0.8)
dat2 <- genFactor(10000, 4, nlevel=4, rho= - 0.1, seed=dat[,6])
dat3 <- genFactor(10000, 4, nlevel=6, rho= -0.2, seed=dat2[,4])
names(dat1) <- sample(LETTERS, length(names(dat1)))
names(dat2) <- sample(letters, length(names(dat2)))
names(dat3) <- sample(paste0(letters,LETTERS), length(names(dat3)))
mdf <- cbind(dat, dat1)
mdf <- cbind(mdf, dat2)
mdf <- cbind(mdf, dat3)

myF <- list(vars = sample(names(mdf), 4))
genFormula(mdf, myF$vars)

misslist <- sample(names(mdf)[c(-1, -22)], 5)
probs <- 0.12
mdf2 <- MAR.df(mdf, vars=misslist, probs=probs)
dimNA(mdf2)

g <- MCARcheck.df(mdf2[, c(-22, -23, -24, -25, -26)])
results <- g
zz <- summary.MCARcheck(g, print=FALSE)

mdf3 <- MCAR.df(mdf, p=probs)

aa <- summary.MCARcheck(MCARcheck.df(mdf3), print=FALSE)

test_that("MAR data is not MCAR", {
  expect_that(aa[[2]], is_less_than(2))
  expect_that(zz[[2]], is_more_than(100))
  expect_that(dimNA(mdf3)[4], approxto(0.12, tol=0.015))
  expect_that(dimNA(mdf2)[4], approxto(0.12, tol=0.015))
  expect_approxto(dimNA(mdf3)[4], 0.12, tol=0.015)
  expect_approxto(dimNA(mdf2)[4], 0.12, tol=0.015)
  expect_approxto(dimNA(mdf2)[4], 0.12, tol=0.00015)
  expect_approxto(dimNA(mdf2)[4], 0.12, tol=1)
})

