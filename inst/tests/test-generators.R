################################################################################
# Test generators
################################################################################

context("Simple Random Correlation Matrices")


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
  expect_is(covmat, "data.frame")
  
})

err <- sum(abs(cor(covmat)[2:ncol(covmat),1] - struc$rho[2:length(struc$rho)])) / length(struc$rho)
tol <- sum(abs(struc$rho[2:length(struc$rho)]) *.05)

test_that("Correlations are correct", {
  expect_that(err, is_less_than(tol))
})



struc2 <- list(dist = c("norm", "binom", "chisq", "pois", "unif", 
                       "weibull", "gamma"), 
              rho = c(-.05, -.4, 0.3, 0.9, .03, -.6, -.2),
              names = c("score", "accept", "score2", "days", "days2", 
                        "luck", "time"))

covmat <- genNumeric(N, pattern=struc2)

test_that("Names get passed properly", {
  expect_equivalent(names(covmat), struc2$names)
  expect_is(covmat, "data.frame")
  
})

context("User specified Seed is Correct")

N <- 1000
K <- 6
RHO1 <- 0.3

S1 <- rnorm(N)
S2 <- runif(N)

struc3 <- list(dist = c("norm", "binom", "chisq", "pois", "unif", 
                        "weibull", "gamma"), 
               rho = c(-.05, -.4, 0.3, 0.9, .03, -.6, -.2),
               names = c("score", "accept", "score2", "days", "days2", 
                         "luck", "time"),
               seed = c(runif(N), rnorm(N), rpois(N, 7), rpois(N, 3), rgamma(N, shape=2), 
                        runif(N)))

covmat <- genNumeric(N, pattern=struc3)

covmat <- genNumeric(N, K, rho=RHO1, seed=S1)
covmat2 <- genNumeric(N, K, rho=RHO1, seed=S2)
test_that("Function appropriately switches between fixed and variable seed",{
  
})

context("Speed")


N <- 98756
K <- 25
P <- 0.43278456
RHO2 <- -0.24
covmat <- genNumeric(N, K, rho=RHO1)

test_that("Function executes in reasonable length of time", {
  expect_that(genNumeric(1000, 25, rho=0.3), takes_less_than(0.5))
  expect_that(genNumeric(500, 40, rho=0.3), takes_less_than(0.5))
  expect_that(genNumeric(87953, 8, rho=0.3), takes_less_than(0.5))
  expect_that(genNumeric(87953, 80, rho=0.3), takes_less_than(1))
})

context("Generate simple correlated factors")

N <- 10000
K <- 30
LEVS <- 12
RHO1 <- 0.2

test <- genFactor(N, K, nlevel=LEVS, rho=RHO1)


test_that("genFactor outputs data in expected size and shape", {
  expect_is(test, "data.frame")
  expect_equal(nrow(test), N)
  expect_equal(ncol(test), K)
  expect_equal(length(table(test[1])), NLEVEL)
})


context("Are the variables in the data.frame properly related?")


chiTest <- function(i, j, data) {chisq.test(data[,i], data[,j])$p.value}
chiP <- Vectorize(chiTest, vectorize.args=list("i", "j"))
Cresults <- outer(1:K, 1:K, chiP, data=test)

gTest <- function(i, j, data) {gammaGK(data[, i], data[, j])$g}
gP <- Vectorize(gTest, vectorize.args=list("i", "j"))
Gresults <- outer(1:K, 1:K, gP, data=test)

gTestSE <- function(i, j, data) {gammaGK(data[, i], data[, j])$se}
gPse <- Vectorize(gTestSE, vectorize.args=list("i", "j"))
GSEresults <- outer(1:K, 1:K, gPse, data=test)

chiSQ <- rep(NA, K-1)

for(i in 2:ncol(Cresults)){
  chiSQ[i-1] <- Cresults[i-1, i] < 0.05
}

summary(chiSQ)

gMAT <- cbind(rep(NA, K-1),rep(NA, K-1)) 

for(i in 2:ncol(Cresults)){
  gMAT[i-1, 1] <- Gresults[i-1, i] + 2* GSEresults[i-1, i]
  gMAT[i-1, 2] <- Gresults[i-1, i] - 2* GSEresults[i-1, i]
}

tol <- RHO1 * .01
test2 <- abs(round(gMAT[,1], digits=3) - RHO1) > tol


test_that("Bivariate relationships exist and magnitude is correct", {
  expect_equivalent(length(chiSQ), table(chiSQ)["TRUE"])
  expect_that(table(test2)["TRUE"][[1]], is_more_than(round(.9*length(test2))))
})

context("Speed tests")

test_that("Function is not slow", {
  expect_that(genFactor(N, K, nlevel=LEVS, rho=RHO1), takes_less_than(1))
})

context("Generate user specified correlated factors")


